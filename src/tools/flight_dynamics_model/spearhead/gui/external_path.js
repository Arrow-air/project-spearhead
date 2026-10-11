import { THREE } from "nicegui-scene";

export default class FlightPath {
  create_mesh(points) {
    const geometry = new THREE.BufferGeometry();
    geometry.setAttribute("position", new THREE.BufferAttribute(new Float32Array(121 * 3), 3));
    this.mesh = new THREE.Line(geometry, new THREE.LineBasicMaterial({
      color: "#22d3ee", transparent: true, opacity: 0.8,
    }));
    this.mesh.frustumCulled = false;
    this.set_points(points);
    return this.mesh;
  }

  set_points(points) {
    const recent = points.slice(-121);
    const position = this.mesh.geometry.attributes.position;
    recent.forEach((point, i) => position.setXYZ(i, ...point));
    position.needsUpdate = true;
    this.mesh.geometry.setDrawRange(0, recent.length);
  }
}
