import { expect, test } from '@playwright/test'

test('loads the real Parkinson atlas and supports core exploration controls', async ({ page }) => {
  await page.goto('./')
  await expect(page.getByRole('heading', { name: /From Genetic Risk/ })).toBeVisible()
  await expect(page.locator('canvas')).toBeVisible()
  await expect(page.getByLabel('Disease', { exact: true })).toHaveValue('parkinson-disease')
  await page.getByLabel('Metric', { exact: true }).selectOption('fdr_p')
  await expect(page.getByLabel('FDR significance legend')).toContainText('0 of 138')
  await page.getByLabel('Metric', { exact: true }).selectOption('spatial_robustness')
  await expect(page.getByLabel('Spatial robustness legend')).toContainText('0 regions meet the joint')
  await page.getByLabel('Metric', { exact: true }).selectOption('validation_score')
  await expect(page.getByLabel('Independent validation legend')).toContainText('Independent validation data')
  await page.getByLabel('Metric', { exact: true }).selectOption('agreement')
  await expect(page.getByLabel('Regional agreement legend')).toContainText('median split')
  await page.getByLabel('Metric', { exact: true }).selectOption('z_score')
  await page.getByLabel('Search brain region', { exact: true }).fill('Putamen L')
  await expect(page.getByLabel('Details for Putamen L')).toBeVisible()
  await expect(page.getByLabel('Details for Putamen L')).toContainText('AAL3v1:77')
  await expect(page.getByLabel('Details for Putamen L')).toContainText('Stage 7 · Spatial sensitivity')
  await expect(page.getByLabel('Details for Putamen L')).toContainText('Stage 8 · Independent validation data')
  await expect(page.getByLabel('Details for Putamen L')).toContainText('Stage 9 · Biological interpretation')
  await expect(page.getByLabel('Details for Putamen L')).toContainText('Top weighted contributors')
  await expect(page.getByLabel('Details for Putamen L')).toContainText('Lput')
  await page.getByRole('button', { name: 'Why is this region highlighted?' }).click()
  await expect(page.getByLabel('Why this region is highlighted')).toContainText('Parkinson genes represented')
  await expect(page.getByLabel('Why this region is highlighted')).toContainText('Evidence:')
  await expect(page.getByLabel('Why this region is highlighted')).toContainText('Interpretation:')
  await page.getByRole('button', { name: 'Show null distribution' }).click()
  await expect(page.getByText(/individual permutation draws were not retained/i)).toBeVisible()
  await page.getByRole('button', { name: 'Reset camera' }).click({ force: true })
  await expect(page.getByRole('link', { name: 'Download regional results' })).toHaveAttribute('download', '')
})

test('atlas validation mode exposes deterministic region IDs', async ({ page }) => {
  await page.goto('./?mode=atlas')
  await expect(page.getByText('ATLAS REGION VALIDATION', { exact: true })).toBeVisible()
  await page.getByRole('button', { name: /Putamen L/ }).click()
  await expect(page.getByLabel('Details for Putamen L')).toContainText('AAL3v1:77')
})

test('WebGL parcel supports direct hover and click hit testing', async ({ page }, testInfo) => {
  test.skip(testInfo.project.name !== 'desktop', 'Mouse hit-testing is covered in desktop mode')
  await page.goto('./?mode=atlas')
  const canvas = page.locator('canvas')
  await expect(canvas).toBeVisible()
  await canvas.scrollIntoViewIfNeeded()
  const box = await canvas.boundingBox()
  expect(box).not.toBeNull()
  let hit = false
  for (const yFraction of [0.3, 0.4, 0.5, 0.6, 0.7]) {
    for (const xFraction of [0.4, 0.5, 0.6, 0.7]) {
      await page.mouse.move(box!.x + box!.width * xFraction, box!.y + box!.height * yFraction)
      if (await page.locator('.tooltip').count()) {
        hit = true
        await page.mouse.down()
        await page.mouse.up()
        break
      }
    }
    if (hit) break
  }
  expect(hit).toBe(true)
  await expect(page.locator('.region-panel')).toBeVisible()
})

test('mobile canvas accepts touch pointer gestures', async ({ page }, testInfo) => {
  test.skip(testInfo.project.name !== 'mobile', 'Touch gestures are covered in mobile mode')
  await page.goto('./')
  const canvas = page.locator('canvas')
  await expect(canvas).toBeVisible()
  await canvas.dispatchEvent('pointerdown', { pointerId: 7, pointerType: 'touch', clientX: 180, clientY: 360, isPrimary: true })
  await canvas.dispatchEvent('pointermove', { pointerId: 7, pointerType: 'touch', clientX: 230, clientY: 390, isPrimary: true })
  await canvas.dispatchEvent('pointerup', { pointerId: 7, pointerType: 'touch', clientX: 230, clientY: 390, isPrimary: true })
  await expect(canvas).toBeVisible()
})
