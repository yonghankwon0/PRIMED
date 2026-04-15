"""
Step 3 (v2): Reader별 Radiomics Feature 추출 — Parallel + CT Windowing
=====================================================================
개선사항:
  1. CT windowing: [-1000, 400] HU clipping (공기/뼈 artifact 제거)
  2. Multiprocessing: 64코어 병렬 처리

사용법: python 03_extract_radiomics_parallel.py [출력경로] [n_workers]
예시:   python 03_extract_radiomics_parallel.py ./output 60
"""

import pylidc as pl
import numpy as np
import pandas as pd
import SimpleITK as sitk
from radiomics import featureextractor
import warnings
import sys
import os
from pathlib import Path
from multiprocessing import Pool, cpu_count
import time

warnings.filterwarnings("ignore")

# =============================================================
# 설정
# =============================================================
OUTPUT_DIR = Path(sys.argv[1]) if len(sys.argv) > 1 else Path("./output")
N_WORKERS = int(sys.argv[2]) if len(sys.argv) > 2 else max(1, cpu_count() - 4)
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

MIN_ANNOTATIONS = 3
CT_MIN = -1000  # 공기 이하 제거
CT_MAX = 400    # 뼈 이상 제거

print(f"=== Parallel Radiomics Extraction ===")
print(f"  Workers: {N_WORKERS}")
print(f"  CT windowing: [{CT_MIN}, {CT_MAX}] HU")
print(f"  Output: {OUTPUT_DIR}")
print()

# =============================================================
# Worker function: 한 scan 처리
# =============================================================
def process_scan(scan_idx):
    """하나의 scan에서 모든 결절의 reader별 radiomics 추출"""
    try:
        # 각 worker에서 독립적으로 extractor 생성
        extractor = featureextractor.RadiomicsFeatureExtractor(
            binWidth=25,
            resampledPixelSpacing=[1, 1, 1],
            interpolator="sitkBSpline",
            enableCExtensions=True,
        )
        extractor.enableAllFeatures()

        # Scan 로드
        scans = pl.query(pl.Scan).all()
        if scan_idx >= len(scans):
            return []
        scan = scans[scan_idx]

        try:
            vol = scan.to_volume().astype(np.float64)
        except Exception:
            return []

        # CT windowing: [-1000, 400] HU clipping
        vol = np.clip(vol, CT_MIN, CT_MAX)

        scan_shape = vol.shape
        spacing = [
            scan.pixel_spacing,
            scan.pixel_spacing,
            scan.slice_thickness or scan.slice_spacing,
        ]

        # SimpleITK image 생성 (clipped volume)
        sitk_image = sitk.GetImageFromArray(vol.T)
        sitk_image.SetSpacing(spacing)

        # Annotation 클러스터링
        try:
            clusters = scan.cluster_annotations(verbose=False)
        except Exception:
            return []

        rows = []
        for cluster in clusters:
            if len(cluster) < MIN_ANNOTATIONS:
                continue

            for reader_idx, ann in enumerate(cluster, start=1):
                try:
                    mask_arr = ann.boolean_mask(pad=0)
                    bbox = ann.bbox()
                    full_mask = np.zeros(scan_shape, dtype=np.uint8)
                    full_mask[bbox] = mask_arr.astype(np.uint8)

                    mask_img = sitk.GetImageFromArray(full_mask.T)
                    mask_img.SetSpacing(spacing)

                    result = extractor.execute(sitk_image, mask_img)
                    features = {
                        k: float(v)
                        for k, v in result.items()
                        if not k.startswith("diagnostics_")
                    }

                    if not features:
                        continue

                    row = {
                        "scan_idx": scan_idx,
                        "patient_id": scan.patient_id,
                        "reader": reader_idx,
                        "n_readers": len(cluster),
                        "subtlety": ann.subtlety,
                        "malignancy": ann.malignancy,
                    }
                    row.update(features)
                    rows.append(row)

                except Exception:
                    continue

        return rows

    except Exception:
        return []


# =============================================================
# Main
# =============================================================
def main():
    t0 = time.time()

    # Scan 수 확인
    scans = pl.query(pl.Scan).all()
    n_scans = len(scans)
    print(f"Total scans: {n_scans}")
    print(f"Starting parallel extraction with {N_WORKERS} workers...")
    print()

    # 병렬 실행
    scan_indices = list(range(n_scans))

    all_rows = []
    completed = 0

    with Pool(processes=N_WORKERS) as pool:
        for result in pool.imap_unordered(process_scan, scan_indices, chunksize=1):
            completed += 1
            all_rows.extend(result)
            if completed % 50 == 0:
                elapsed = time.time() - t0
                rate = completed / elapsed * 60
                eta = (n_scans - completed) / (rate / 60) if rate > 0 else 0
                print(
                    f"  [{completed}/{n_scans}] rows: {len(all_rows)}, "
                    f"rate: {rate:.0f} scans/min, ETA: {eta/60:.1f}min"
                )

    elapsed = time.time() - t0
    print(f"\nExtraction done in {elapsed/60:.1f} min")

    # =============================================================
    # 결절 ID 할당 (같은 scan의 같은 cluster → 같은 nodule_id)
    # =============================================================
    df = pd.DataFrame(all_rows)
    if len(df) == 0:
        print("ERROR: No data extracted!")
        return

    # scan_idx + reader 그룹 기반으로 nodule_id 부여
    # 같은 scan에서 같은 reader 수를 가진 연속 그룹이 같은 결절
    df = df.sort_values(["scan_idx", "n_readers", "reader"]).reset_index(drop=True)

    nodule_id = 0
    prev_scan = -1
    prev_nreaders = -1
    reader_count = 0
    nodule_ids = []

    for _, row in df.iterrows():
        if row["scan_idx"] != prev_scan:
            nodule_id += 1
            reader_count = 1
            prev_scan = row["scan_idx"]
            prev_nreaders = row["n_readers"]
        elif reader_count >= prev_nreaders:
            nodule_id += 1
            reader_count = 1
            prev_nreaders = row["n_readers"]
        else:
            reader_count += 1
        nodule_ids.append(nodule_id)

    df["nodule_id"] = nodule_ids

    # =============================================================
    # 저장: per-reader CSV
    # =============================================================
    per_reader_path = OUTPUT_DIR / "radiomics_per_reader.csv"
    df.to_csv(per_reader_path, index=False)
    print(f"\n=== Per-reader CSV: {per_reader_path} ===")
    print(f"  Rows: {len(df)}, Nodules: {df['nodule_id'].nunique()}")
    feature_cols_list = [c for c in df.columns if c.startswith("original_")]
    print(f"  Features: {len(feature_cols_list)}")

    # =============================================================
    # Consensus-Disagreement 계산
    # =============================================================
    meta_cols = [
        "scan_idx", "nodule_id", "patient_id", "reader", "n_readers",
        "subtlety", "malignancy",
    ]
    feature_cols = [c for c in df.columns if c not in meta_cols]

    cd_rows = []
    for nid, grp in df.groupby("nodule_id"):
        row = {
            "nodule_id": nid,
            "patient_id": grp["patient_id"].iloc[0],
            "n_readers": grp["n_readers"].iloc[0],
            "malignancy_mean": grp["malignancy"].mean(),
        }
        for feat in feature_cols:
            vals = grp[feat].dropna()
            if len(vals) >= 2:
                row[f"{feat}_cons"] = vals.mean()
                row[f"{feat}_disp"] = vals.std()
            else:
                row[f"{feat}_cons"] = vals.mean() if len(vals) > 0 else np.nan
                row[f"{feat}_disp"] = 0.0
        cd_rows.append(row)

    cd_df = pd.DataFrame(cd_rows)
    cd_path = OUTPUT_DIR / "radiomics_cd.csv"
    cd_df.to_csv(cd_path, index=False)
    print(f"\n=== Consensus-Disagreement CSV: {cd_path} ===")
    print(f"  Nodules: {len(cd_df)}")
    cons_count = len([c for c in cd_df.columns if c.endswith("_cons")])
    disp_count = len([c for c in cd_df.columns if c.endswith("_disp")])
    print(f"  Consensus features: {cons_count}")
    print(f"  Disagreement features: {disp_count}")

    total_time = time.time() - t0
    print(f"\n=== Total time: {total_time/60:.1f} min ===")


if __name__ == "__main__":
    main()
