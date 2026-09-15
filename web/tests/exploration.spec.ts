import { expect, test } from '@playwright/test'

test('loads the real multi-disease atlas and supports metric controls', async ({ page }, testInfo) => {
  test.setTimeout(120_000)
  await page.goto('./')
  await expect(page.getByRole('heading', { level: 1, name: /From Genetic Risk to Spatial Brain Vulnerability/ })).toBeVisible()
  await expect(page.locator('canvas').first()).toBeVisible({ timeout: 30_000 })
  await expect(page.getByLabel('Disease', { exact: true })).toHaveValue('parkinson')
  if (testInfo.project.name === 'mobile') {
    await page.getByLabel('Disease', { exact: true }).selectOption('alzheimer')
    await expect(page.getByText(/ALZHEIMER DISEASE · AAL3 · 138 REGIONS/)).toBeVisible()
    return
  }
  await page.getByLabel('Metric', { exact: true }).selectOption('fdr_p')
  await expect(page.getByLabel('FDR significance legend')).toContainText('0 of 138')
  await page.getByLabel('Metric', { exact: true }).selectOption('spatial_robustness')
  await expect(page.getByLabel('Spatial robustness legend')).toContainText('0 regions meet the joint')
  await page.getByLabel('Metric', { exact: true }).selectOption('validation_score')
  await expect(page.getByLabel('Independent validation legend')).toContainText('Independent validation data')
  await page.getByLabel('Metric', { exact: true }).selectOption('agreement')
  await expect(page.getByLabel('Regional agreement legend')).toContainText('median split')
  await page.getByLabel('Metric', { exact: true }).selectOption('z_score')
})

test('shows regional evidence, null summary, and downloads', async ({ page }, testInfo) => {
  test.skip(testInfo.project.name === 'mobile', 'Detailed disclosures are covered on desktop')
  test.setTimeout(120_000)
  await page.goto('./')
  await expect(page.getByLabel('Search brain region', { exact: true })).toBeVisible()
  await page.getByLabel('Search brain region', { exact: true }).fill('Putamen L')
  await expect(page.getByLabel('Details for Putamen L')).toBeVisible()
  await expect(page.getByLabel('Details for Putamen L')).toContainText('AAL3v1:77')
  await expect(page.getByLabel('Details for Putamen L')).toContainText('Stage 7 · Spatial sensitivity')
  await expect(page.getByLabel('Details for Putamen L')).toContainText('Independent validation data')
  await expect(page.getByLabel('Details for Putamen L')).toContainText('Biological interpretation')
  await expect(page.getByLabel('Details for Putamen L')).toContainText('Lput')
  await page.getByLabel('Search brain region', { exact: true }).fill('OFCant R')
  await expect(page.getByLabel('Details for OFCant R')).toContainText('Top weighted contributors')
  await page.getByRole('button', { name: 'Why is this region highlighted?' }).click()
  await expect(page.getByLabel('Why this region is highlighted')).toContainText('Disease genes represented')
  await expect(page.getByLabel('Why this region is highlighted')).toContainText('Evidence:')
  await expect(page.getByLabel('Why this region is highlighted')).toContainText('Interpretation:')
  await page.getByRole('button', { name: 'Show null distribution' }).click()
  await expect(page.getByText(/individual permutation draws were not retained/i)).toBeVisible()
  await page.getByRole('button', { name: 'Reset camera' }).click({ force: true })
  await expect(page.getByRole('link', { name: 'Download regional results' })).toHaveAttribute('download', '')
})

test('switches disease and shows comparison', async ({ page }, testInfo) => {
  test.skip(testInfo.project.name === 'mobile', 'Mobile disease switching is covered in the metric-controls test')
  await page.goto('./')
  await expect(page.getByLabel('Disease', { exact: true })).toHaveValue('parkinson')
  await page.getByLabel('Disease', { exact: true }).selectOption('alzheimer')
  await expect(page.getByText(/ALZHEIMER DISEASE · AAL3 · 138 REGIONS/)).toBeVisible()
  await expect(page.getByRole('heading', { name: 'Compare Diseases' })).toBeVisible()
})

test('atlas validation mode exposes deterministic region IDs', async ({ page }) => {
  await page.goto('./?mode=atlas')
  await expect(page.getByText('ATLAS REGION VALIDATION', { exact: true })).toBeVisible()
  await page.getByRole('button', { name: /Putamen L/ }).click()
  await expect(page.getByLabel('Details for Putamen L')).toContainText('AAL3v1:77')
})

test('mobile canvas accepts touch pointer gestures', async ({ page }, testInfo) => {
  test.skip(testInfo.project.name !== 'mobile', 'Touch gestures are covered in mobile mode')
  await page.goto('./')
  const canvas = page.locator('canvas').first()
  await expect(canvas).toBeVisible({ timeout: 30_000 })
  await canvas.dispatchEvent('pointerdown', { pointerId: 7, pointerType: 'touch', clientX: 180, clientY: 360, isPrimary: true })
  await canvas.dispatchEvent('pointermove', { pointerId: 7, pointerType: 'touch', clientX: 230, clientY: 390, isPrimary: true })
  await canvas.dispatchEvent('pointerup', { pointerId: 7, pointerType: 'touch', clientX: 230, clientY: 390, isPrimary: true })
  await expect(canvas).toBeVisible()
})
