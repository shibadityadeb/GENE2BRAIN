export function formatValue(value: number | null | undefined, digits = 3) {
  return Number.isFinite(value) ? Number(value).toFixed(digits) : 'N/A'
}

export function formatCount(value: number | null | undefined) {
  return Number.isFinite(value) ? Number(value).toLocaleString('en-US') : 'N/A'
}
