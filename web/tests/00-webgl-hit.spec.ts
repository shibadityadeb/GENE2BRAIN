import { expect, test } from '@playwright/test'

test('WebGL parcel supports direct hover and click hit testing', async ({ page }, testInfo) => {
  test.skip(testInfo.project.name !== 'desktop', 'Mouse hit-testing is covered in desktop mode')
  await page.setViewportSize({ width: 1440, height: 1100 })
  await page.goto('./?mode=atlas')
  const canvas = page.locator('canvas').first()
  await expect(canvas).toBeVisible()
  await canvas.scrollIntoViewIfNeeded()
  await page.waitForTimeout(750)
  const box = await canvas.boundingBox()
  expect(box).not.toBeNull()
  // Validation mode has a deterministic camera; the centre point contains
  // parcel geometry and therefore exercises React Three Fiber ray-casting.
  await page.mouse.move(box!.x + box!.width * 0.55, box!.y + box!.height * 0.5)
  await expect(page.locator('.tooltip')).toBeVisible()
  await page.mouse.down()
  await page.mouse.up()
  await expect(page.locator('.region-panel')).toBeVisible()
})
