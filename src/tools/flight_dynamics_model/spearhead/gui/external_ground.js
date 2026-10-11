import { THREE } from "nicegui-scene";

export default class Ground {
  create_mesh() {
    const group = new THREE.Group();
    const plane = new THREE.Mesh(
      new THREE.PlaneGeometry(6000, 6000),
      new THREE.MeshPhongMaterial({ color: "#3e6259", side: THREE.DoubleSide }),
    );
    plane.position.z = -0.02;
    group.add(plane);
    const grid = new THREE.GridHelper(6000, 240, "#a9c9bd", "#a9c9bd");
    grid.rotateX(Math.PI / 2);
    grid.material.transparent = true;
    grid.material.opacity = 0.22;
    group.add(grid);
    return group;
  }
}
