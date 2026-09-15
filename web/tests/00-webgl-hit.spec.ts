import { expect, test } from '@playwright/test'

test('WebGL parcel supports direct hover and click hit testing', async ({ page }, testInfo) => {
  test.skip(testInfo.project.name !== 'desktop', 'Mouse hit-testing is covered in desktop mode')
  await page.addInitScript(() => window.localStorage.setItem('gene2brain-guide-seen', '1'))
  await page.setViewportSize({ width: 1440, height: 1100 })
  await page.goto('./?mode=atlas')
  const canvas = page.locator('canvas').first()
  await expect(canvas).toBeVisible({ timeout: 30_000 })
  await canvas.scrollIntoViewIfNeeded()
  await page.waitForTimeout(750)
  const box = await canvas.boundingBox()
  expect(box).not.toBeNull()
  // Scan the deterministic atlas viewport for a rendered parcel instead of
  // assuming a single pixel remains covered after geometry refinements.
  const xFractions = [0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8]
  const yFractions = [0.25, 0.35, 0.45, 0.55, 0.65, 0.75]
  let hitPoint: { x: number; y: number } | undefined
  for (const y of yFractions) {
    for (const x of xFractions) {
      const point = { x: box!.x + box!.width * x, y: box!.y + box!.height * y }
      await page.mouse.move(point.x, point.y)
      if (await page.locator('.tooltip').isVisible()) {
        hitPoint = point
        break
      }
    }
    if (hitPoint) break
  }
  expect(hitPoint).toBeDefined()
  await page.mouse.click(hitPoint!.x, hitPoint!.y)
  await expect(page.locator('.region-panel')).toBeVisible()
})
