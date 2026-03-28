# Milo Differential Abundance Analysis — CRC Single-Cell Data

대장암(CRC) 단일세포 RNA-seq 데이터에 Milo 기반 차등 풍부도(DA) 분석을 적용하여,
MSI(현미부수체 불안정성) vs MSS(현미부수체 안정성) 종양 미세환경을 비교합니다.

이 레포지토리는 더 넓은 ARPA-H 관련 CRC 단일세포 분석 프로젝트에서 제가 담당한 MILO 기반 차등 풍부도 분석 모듈을 포함합니다.

---

## Milo란?

[Milo](https://github.com/emdann/milopy)는 단일세포 데이터에 대한 차등 풍부도 검정을 위한 통계적 방법입니다.
사전 정의된 클러스터를 비교하는 대신, k-최근접 이웃(KNN) 이웃 구조에서의 풍부도 변화를 검정함으로써
임의적인 클러스터링 결정에 대한 민감도를 줄입니다.

> Dann et al. (2022) *Differential abundance testing on single-cell data using
> k-nearest neighbor graphs.* Nature Biotechnology.

---

## 스크립트 동작

1. `.h5ad` 단일세포 데이터셋을 불러오고, CSV로부터 컨센서스 클러스터 레이블을 매핑합니다
2. 현미부수체 상태(MSI / MSS)로 세포를 필터링합니다
3. PCA, KNN 그래프 구성, UMAP을 실행합니다
4. Milo DA 파이프라인을 실행합니다 (`make_nhoods` → `count_nhoods` → `DA_nhoods`)
5. logFC, p-value, spatial FDR 기준으로 유의미한 이웃을 필터링합니다
6. MSI 및 MSS 농축 세포 목록을 CSV 파일로 내보냅니다
7. DA 이웃 그래프와 volcano plot을 저장합니다

---

## 사용법

```bash
pip install -r requirements.txt

python MILO.py \
    --input_h5ad  /path/to/B_cell.h5ad \
    --input_csv   /path/to/consensus_clustering_result.csv \
    --output_dir  ./output
```

### 주요 인수

| 인수 | 기본값 | 설명 |
|------|--------|------|
| `--input_h5ad` | 필수 | `.h5ad` 입력 파일 경로 |
| `--input_csv` | 필수 | 컨센서스 클러스터링 CSV 경로 |
| `--output_dir` | `./output` | 출력 파일 디렉토리 |
| `--sample_col` | `sample_id` | 생물학적 반복 실험을 위한 obs 열 |
| `--status_col` | `microsatellite_status` | MSI/MSS 조건을 위한 obs 열 |
| `--alpha` | `0.01` | Spatial FDR 임계값 |
| `--logfc_thr` | `0.2` | logFC 임계값 |
| `--pval_thr` | `0.05` | p-value 임계값 |
| `--seed` | `42` | 랜덤 시드 |

