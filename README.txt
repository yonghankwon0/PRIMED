PRIMED: Joint Loss Estimation with Selective Penalization
for Multi-Reader Imaging Data with Binary Outcomes
=========================================================

폴더 구조
---------

paper/                    논문 컴파일에 필요한 모든 파일
  primed_paper.tex          영문 논문 (메인)
  primed_paper_korean.tex   한국어 번역본
  primed_paper.pdf          컴파일된 영문 PDF
  primed_paper_korean.pdf   컴파일된 한국어 PDF
  sn-jnl.cls               Springer Nature 저널 클래스
  sn-bibliography.bib      참고문헌 데이터베이스
  sn-vancouver-num.bst     참고문헌 스타일 (Vancouver numbered)
  sn-basic.bst             참고문헌 스타일 (basic)
  primed_paper.bbl         컴파일된 참고문헌
  fig_path_monotonicity.pdf   Figure: 정규화 경로 단조성
  fig_scatter_cd.pdf          Figure: 합의-불일치 산점도 (직장암)

figures/                  논문 및 분석에 사용된 그림 파일
  fig_path_monotonicity.pdf   정규화 경로 단조성 (LIDC-IDRI)
  fig_scatter_cd.pdf          합의-불일치 산점도 (직장암, 10개 순서 특징)
  lidc_cd_scatter.pdf         LIDC-IDRI 합의-불일치 산점도
  lidc_cd_top12_scatter.pdf   LIDC-IDRI 상위 12개 특징 산점도
  lidc_cd_all107_summary.pdf  LIDC-IDRI 107개 특징 요약

data/                     입력 데이터
  data.xlsx                 직장암 multi-reader MRI 데이터 (n=149, m=5, K=14)
  lidc_annotations.csv      LIDC-IDRI 폐 결절 주석 데이터 (n=1,390, m=3-4)

code/                     분석 코드
  rectal/                   직장암 데이터 분석
    01_cd_group_lasso.R       Group LASSO 분석 (CD GL, C LASSO)
    02_primed_linear.R        PRIMED 적용 (선형 산포)
    03_primed_quadratic.R     PRIMED 적용 (이차 산포)
    04_bootstrap_stability.R  부트스트랩 안정성 분석
    05_scatter_plot.R         합의-불일치 산점도 생성

  lidc/                     LIDC-IDRI 폐 결절 분석
    01_extract_radiomics.py   PyRadiomics 특징 추출 (60-worker 병렬)
    02_apply_primed.R         PRIMED 적용
    03_cv10_evaluation.R      10-fold nested CV 평가
    04_bootstrap_stability.R  부트스트랩 안정성 분석

  simulation/               시뮬레이션 연구
    01_primed_core.jl         PRIMED Julia 구현 (K<=50)
    01_primed_core_highdim.jl PRIMED Julia 구현 (K=100)
    02_run_single_job.R       단일 시뮬레이션 작업 (R wrapper)
    03_run_parallel.sh        병렬 실행 스크립트
    04_aggregate_results.R    결과 집계 스크립트

results/                  분석 결과
  simulation/               시뮬레이션 결과
    sim_all_with_sd_results.csv   최종 집계 결과 (40 시나리오 x 3 방법)
    sim_all_selfreq.csv           변수 선택 빈도
    result_full_job01-40.csv      개별 작업별 결과
    selfreq_full_job01-40.csv     개별 작업별 선택 빈도

  lidc/                     LIDC-IDRI 결과
    radiomics_cd.csv              합의-불일치 radiomics 행렬
    radiomics_per_reader.csv      판독자별 radiomics 행렬
    cv10v2_auc_results.csv        10-fold CV AUC 결과
    cv10v2_selection_frequency.csv  변수 선택 빈도
    cv10v2_full_results.rds       전체 CV 결과 (R 객체)
    primed_radiomics_results.rds  PRIMED 적합 결과 (R 객체)


논문 컴파일
-----------
cd paper/
tectonic primed_paper.tex
tectonic primed_paper_korean.tex


핵심 방법
---------
- 결합 손실: L_disp/L_disp^0 + L_out/L_out^0
- 잔차 불일치: S_tilde = S - alpha_0 - alpha * X_bar
- Group LASSO: (beta_k, gamma_k)에만 벌칙, alpha_k는 비벌칙
- 데이터 기반 lambda 경로: lambda_max로부터 로그 등간격
