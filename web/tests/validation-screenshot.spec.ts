import { expect, test } from '@playwright/test'

test('capture AAL3 region validation view', async ({ page }, testInfo) => {
  test.skip(testInfo.project.name !== 'desktop', 'Canonical screenshot uses the desktop viewport')
  await page.setViewportSize({ width: 1440, height: 1100 })
  await page.goto('./?mode=atlas')
  await expect(page.getByText('ATLAS REGION VALIDATION', { exact: true })).toBeVisible()
  await expect(page.locator('canvas').first()).toBeVisible()
  await page.waitForTimeout(1500)
  await page.locator('.brain-stage').screenshot({ path: '../results/figures/web_atlas_region_validation.png' })
})
