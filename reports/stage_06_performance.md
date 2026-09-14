# Stage 6 performance

- Regions: 138
- AHBA genes: 15,632
- Eligible non-Parkinson background genes: 15,509
- Represented genes: broad=123, stringent=34, weighted=123
- Permutations per analysis: 10,000
- Total runtime: 7.62 seconds
- Approximate peak resident memory: 789.6 MB
- Batch size: 250 permutations
- Memory strategy: random indices are stored as compact int32 matrices; null scores
  are stored as three two-dimensional permutation × region float32 arrays. Temporary
  region × batch × gene arrays are released each batch; no full three-dimensional
  permutation object is retained.
- Python: 3.12.7
- Platform: macOS-26.1-arm64-arm-64bit
- Processor: arm
- Logical CPU count: 8
- Random seed: 20260913
