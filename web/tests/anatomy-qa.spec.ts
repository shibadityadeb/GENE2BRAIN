import { expect, test } from '@playwright/test'

test('capture anatomy-only quality-control view', async ({ page }, testInfo) => {
  test.skip(testInfo.project.name !== 'desktop', 'Canonical anatomy screenshot uses desktop viewport')
  await page.setViewportSize({ width: 1440, height: 1100 })
  await page.goto('./')
  await page.getByLabel('Metric', { exact: true }).selectOption('anatomy')
  await expect(page.getByLabel('Anatomy only legend')).toBeVisible()
  await page.waitForTimeout(1200)
  await page.locator('.brain-stage').screenshot({ path: '../results/figures/web_brain_anatomical.png' })
})

test('capture real-data brain and atlas mapping views', async ({ page }, testInfo) => {
  test.skip(testInfo.project.name !== 'desktop', 'Canonical screenshots use desktop viewport')
  await page.setViewportSize({ width: 1440, height: 1100 })
  await page.goto('./')
  await expect(page.locator('.brain-stage').first()).toBeVisible()
  await page.waitForTimeout(1200)
  await page.locator('.brain-stage').first().screenshot({ path: '../results/figures/web_brain_gene2brain.png' })
  await page.goto('./?mode=atlas')
  await expect(page.getByLabel('Atlas validation legend')).toBeVisible()
  await page.waitForTimeout(1200)
  await page.locator('.brain-stage').first().screenshot({ path: '../results/figures/web_atlas_region_validation.png' })
})

test('check left, right and superior anatomical orientation', async ({ page }, testInfo) => {
  test.skip(testInfo.project.name !== 'desktop', 'Anatomical orientation screenshots use desktop viewport')
  await page.setViewportSize({ width: 1440, height: 1100 })
  await page.goto('/')
  await page.getByLabel('Metric', { exact: true }).selectOption('anatomy')
  for (const view of ['Left', 'Right', 'Superior']) {
    await page.getByRole('button', { name: view, exact: true }).last().click()
    await page.waitForTimeout(950)
    await page.locator('.brain-stage').screenshot({ path: `../results/figures/web_brain_${view.toLowerCase()}_anatomical.png` })
  }
})
