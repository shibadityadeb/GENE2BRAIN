import { Canvas, type ThreeEvent, useFrame, useThree } from '@react-three/fiber'
import { useEffect, useMemo, useRef } from 'react'
import {
  BufferAttribute, BufferGeometry, Color, DoubleSide, Mesh, MeshBasicMaterial,
  MeshStandardMaterial, TOUCH, Vector3,
} from 'three'
import { OrbitControls as OrbitControlsImpl } from 'three/addons/controls/OrbitControls.js'
import { colorForRegion } from '../colors'
import type { AtlasGeometry, Hemisphere, Metric, RegionRecord, ViewPreset } from '../types'

interface Props {
  geometryData: AtlasGeometry
  regions: RegionRecord[]
  metric: Metric
  hemisphere: Hemisphere
  selected: RegionRecord | null
  focusNonce: number
  viewPreset: ViewPreset
  validationMode: boolean
  zDomain: [number, number]
  observedDomain: [number, number]
  fdrThreshold: number
  onHover: (region: RegionRecord | null, point?: { x: number; y: number }) => void
  onSelect: (region: RegionRecord) => void
  onReset: () => void
}

interface CombinedGeometry {
  geometry: BufferGeometry
  faceRegions: number[]
  regionGeometry: Map<number, BufferGeometry>
  regionVertexRanges: Map<number, [number, number]>
  center: Vector3
}

function combineGeometry(
  source: AtlasGeometry,
  allowed: Set<number>,
): CombinedGeometry {
  const positions: number[] = []
  const colors: number[] = []
  const indices: number[] = []
  const faceRegions: number[] = []
  const regionGeometry = new Map<number, BufferGeometry>()
  const regionVertexRanges = new Map<number, [number, number]>()
  const boundsGeometry = new BufferGeometry()

  source.regions.forEach((parcel) => {
    if (!allowed.has(parcel.region_id)) return
    const offset = positions.length / 3
    positions.push(...parcel.positions)
    const vertexCount = parcel.positions.length / 3
    regionVertexRanges.set(parcel.region_id, [offset, vertexCount])
    for (let index = 0; index < vertexCount; index += 1) colors.push(1, 1, 1)
    indices.push(...parcel.indices.map((index) => index + offset))
    for (let index = 0; index < parcel.indices.length / 3; index += 1) faceRegions.push(parcel.region_id)
  })
  boundsGeometry.setAttribute('position', new BufferAttribute(new Float32Array(positions), 3))
  boundsGeometry.computeBoundingBox()
  const center = new Vector3()
  boundsGeometry.boundingBox!.getCenter(center)
  for (let index = 0; index < positions.length; index += 3) {
    positions[index] -= center.x
    positions[index + 1] -= center.y
    positions[index + 2] -= center.z
  }
  const geometry = new BufferGeometry()
  geometry.setAttribute('position', new BufferAttribute(new Float32Array(positions), 3))
  geometry.setAttribute('color', new BufferAttribute(new Float32Array(colors), 3))
  geometry.setIndex(indices)
  geometry.computeVertexNormals()
  geometry.computeBoundingSphere()

  source.regions.forEach((parcel) => {
    if (!allowed.has(parcel.region_id)) return
    const parcelPositions = [...parcel.positions]
    for (let index = 0; index < parcelPositions.length; index += 3) {
      parcelPositions[index] -= center.x
      parcelPositions[index + 1] -= center.y
      parcelPositions[index + 2] -= center.z
    }
    const parcelGeometry = new BufferGeometry()
    parcelGeometry.setAttribute('position', new BufferAttribute(new Float32Array(parcelPositions), 3))
    parcelGeometry.setIndex(parcel.indices)
    parcelGeometry.computeVertexNormals()
    regionGeometry.set(parcel.region_id, parcelGeometry)
  })
  boundsGeometry.dispose()
  return { geometry, faceRegions, regionGeometry, regionVertexRanges, center }
}

const PRESET_POSITIONS: Record<ViewPreset, Vector3> = {
  reset: new Vector3(150, 80, 180),
  anterior: new Vector3(0, 10, 240),
  posterior: new Vector3(0, 10, -240),
  superior: new Vector3(0, 240, 0.01),
  inferior: new Vector3(0, -240, 0.01),
}

function CameraRig({ selected, center, focusNonce, viewPreset }: {
  selected: RegionRecord | null
  center: Vector3
  focusNonce: number
  viewPreset: ViewPreset
}) {
  const { camera, gl } = useThree()
  const controls = useMemo(() => {
    const instance = new OrbitControlsImpl(camera, gl.domElement)
    instance.enableDamping = true
    instance.dampingFactor = 0.08
    instance.minDistance = 60
    instance.maxDistance = 400
    instance.screenSpacePanning = true
    instance.touches.ONE = TOUCH.ROTATE
    instance.touches.TWO = TOUCH.DOLLY_PAN
    return instance
  }, [camera, gl.domElement])
  const destinationPosition = useRef(PRESET_POSITIONS.reset.clone())
  const destinationTarget = useRef(new Vector3())
  const animating = useRef(true)

  useEffect(() => {
    if (selected && focusNonce > 0) {
      const [x, y, z] = selected.centroid_mni
      destinationTarget.current.set(x - center.x, z - center.y, y - center.z)
      const direction = camera.position.clone().sub(controls.target).normalize()
      destinationPosition.current.copy(destinationTarget.current).add(direction.multiplyScalar(115))
    } else {
      destinationTarget.current.set(0, 0, 0)
      destinationPosition.current.copy(PRESET_POSITIONS[viewPreset])
    }
    animating.current = true
  }, [camera, center, focusNonce, selected, viewPreset])

  useFrame(() => {
    controls.update()
    if (!animating.current) return
    camera.position.lerp(destinationPosition.current, 0.09)
    controls.target.lerp(destinationTarget.current, 0.09)
    if (camera.position.distanceTo(destinationPosition.current) < 0.2 && controls.target.distanceTo(destinationTarget.current) < 0.2) {
      animating.current = false
    }
  })
  useEffect(() => () => controls.dispose(), [controls])

  return <primitive object={controls} />
}

function AtlasMesh(props: Props) {
  const records = useMemo(() => new Map(props.regions.map((region) => [region.region_id, region])), [props.regions])
  const allowed = useMemo(() => new Set(props.regions
    .filter((region) => props.hemisphere === 'whole' || region.hemisphere === props.hemisphere || region.hemisphere === 'B')
    .map((region) => region.region_id)), [props.hemisphere, props.regions])
  const combined = useMemo(() => combineGeometry(props.geometryData, allowed), [allowed, props.geometryData])
  const material = useMemo(() => new MeshStandardMaterial({
    vertexColors: true, roughness: 0.72, metalness: 0.04, side: DoubleSide,
    transparent: props.validationMode, opacity: props.validationMode ? 0.9 : 1,
  }), [props.validationMode])
  const highlightMaterial = useMemo(() => new MeshBasicMaterial({
    color: new Color('#ffffff'), transparent: true, opacity: 0.32, depthTest: false, side: DoubleSide,
  }), [])
  const selectedMaterial = useMemo(() => new MeshBasicMaterial({
    color: new Color('#ffe08a'), transparent: true, opacity: 0.72, depthTest: false,
    side: DoubleSide, wireframe: true,
  }), [])
  const hoverId = useRef<number | null>(null)
  const hoverMesh = useRef<Mesh>(null)

  useEffect(() => () => {
    combined.geometry.dispose()
    combined.regionGeometry.forEach((geometry) => geometry.dispose())
  }, [combined])
  useEffect(() => {
    const colorAttribute = combined.geometry.getAttribute('color') as BufferAttribute
    combined.regionVertexRanges.forEach(([offset, count], regionId) => {
      const color = colorForRegion(
        records.get(regionId)!, props.metric,
        { zDomain: props.zDomain, observedDomain: props.observedDomain, fdrThreshold: props.fdrThreshold },
        props.validationMode,
      )
      for (let index = offset; index < offset + count; index += 1) colorAttribute.setXYZ(index, color.r, color.g, color.b)
    })
    colorAttribute.needsUpdate = true
  }, [combined, props.fdrThreshold, props.metric, props.observedDomain, props.validationMode, props.zDomain, records])
  useEffect(() => () => { material.dispose(); highlightMaterial.dispose(); selectedMaterial.dispose() }, [highlightMaterial, material, selectedMaterial])

  const setHoverGeometry = (regionId: number | null) => {
    hoverId.current = regionId
    if (hoverMesh.current) {
      hoverMesh.current.visible = regionId !== null
      if (regionId !== null) hoverMesh.current.geometry = combined.regionGeometry.get(regionId)!
    }
  }

  const onMove = (event: ThreeEvent<PointerEvent>) => {
    event.stopPropagation()
    if (event.faceIndex == null) return
    const regionId = combined.faceRegions[event.faceIndex]
    if (hoverId.current !== regionId) setHoverGeometry(regionId)
    props.onHover(records.get(regionId) ?? null, { x: event.nativeEvent.clientX, y: event.nativeEvent.clientY })
  }
  const onOut = () => { setHoverGeometry(null); props.onHover(null) }
  const onClick = (event: ThreeEvent<MouseEvent>) => {
    event.stopPropagation()
    if (event.faceIndex == null) return
    const record = records.get(combined.faceRegions[event.faceIndex])
    if (record) props.onSelect(record)
  }

  return (
    <>
      <mesh
        geometry={combined.geometry}
        material={material}
        onPointerMove={onMove}
        onPointerOut={onOut}
        onClick={onClick}
      />
      <mesh ref={hoverMesh} material={highlightMaterial} renderOrder={5} visible={false} />
      {props.selected && allowed.has(props.selected.region_id) && (
        <mesh geometry={combined.regionGeometry.get(props.selected.region_id)} material={selectedMaterial} renderOrder={6} />
      )}
      <CameraRig
        selected={props.selected}
        center={combined.center}
        focusNonce={props.focusNonce}
        viewPreset={props.viewPreset}
      />
    </>
  )
}

export function BrainScene(props: Props) {
  return (
    <Canvas
      aria-label="Interactive AAL3 brain atlas"
      camera={{ position: [150, 80, 180], fov: 42, near: 0.1, far: 1000 }}
      dpr={[1, 1.7]}
      gl={{ antialias: true, powerPreference: 'high-performance' }}
      onDoubleClick={props.onReset}
    >
      <color attach="background" args={['#071013']} />
      <ambientLight intensity={1.7} />
      <directionalLight position={[100, 150, 200]} intensity={2.2} />
      <directionalLight position={[-120, -50, -100]} intensity={0.8} color="#8ac6c9" />
      <AtlasMesh {...props} />
    </Canvas>
  )
}
