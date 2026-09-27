import * as THREE from 'three';
import {OrbitControls} from 'three/addons/controls/OrbitControls.js';
import {GLTFLoader} from 'three/addons/loaders/GLTFLoader.js';
import {mergeGeometries} from 'three/addons/utils/BufferGeometryUtils.js';
const stage=document.querySelector('#stage');
const renderer=new THREE.WebGLRenderer({antialias:true,alpha:false});
renderer.setPixelRatio(Math.min(devicePixelRatio,2));renderer.setClearColor(0x10191f);
renderer.outputColorSpace=THREE.SRGBColorSpace;renderer.toneMapping=THREE.ACESFilmicToneMapping;renderer.toneMappingExposure=1.1;
stage.appendChild(renderer.domElement);renderer.domElement.setAttribute('aria-label','Spearhead 3D model. Drag to rotate, scroll to zoom.');
const scene=new THREE.Scene();
scene.add(new THREE.HemisphereLight(0xe1efff,0x45504a,2.7));
for(const [pos,color,intensity] of [[[2,5,3],0xffffff,3],[[-3,2,-4],0xb9d5ff,2],[[1,-2,1],0xffe4bf,.7]]){const light=new THREE.DirectionalLight(color,intensity);light.position.set(...pos);scene.add(light);}
const camera=new THREE.PerspectiveCamera(35,1,.01,100);
const controls=new OrbitControls(camera,renderer.domElement);controls.enableDamping=true;controls.dampingFactor=.08;controls.minDistance=.05;controls.maxDistance=30;
const names={fuselage:'Fuselage',main_wing:'Main wing · source',tail:'V-tail',motor_mounts_and_booms:'Booms & motor mounts',landing_gear:'Front landing gear',inferred_starboard_outer_wing:'Mirrored outer wing',restored_hidden_port_aft_motor_mount:'Recovered fourth mount',restored_center_wing_ribs:'Center-wing ribs',restored_rear_landing_gear:'Rear landing gear'};
const colors={fuselage:'#cead76',main_wing:'#d8c49b',tail:'#cfb27e',motor_mounts_and_booms:'#657986',landing_gear:'#ad987b',inferred_starboard_outer_wing:'#87b4cf',restored_hidden_port_aft_motor_mount:'#75c3ba',restored_center_wing_ribs:'#75c3ba',restored_rear_landing_gear:'#75c3ba'};
let root,meta,groups=[],meshes=[],selected,view='iso';const center=new THREE.Vector3();let radius=2;
let wire=false;
function clearSelection(){if(selected){selected.material.emissive.setHex(0);selected=null;}document.querySelector('#selection strong').textContent='Click a part to identify it';}
function updateCount(){document.querySelector('#count').textContent=groups.reduce((sum,g)=>sum+(g.object.visible?g.solids:0),0);for(const g of groups)g.checkbox.checked=g.object.visible;clearSelection();}
function all(){for(const g of groups)g.object.visible=true;updateCount();}
function fit(next=view){if(!root)return;view=next;document.querySelectorAll('[data-view]').forEach(b=>b.classList.toggle('active',b.dataset.view===view));
const box=new THREE.Box3();for(const g of groups)if(g.object.visible)box.union(new THREE.Box3().setFromObject(g.object));if(box.isEmpty())box.setFromObject(root);
box.getCenter(center);const size=box.getSize(new THREE.Vector3());radius=size.length()/2;
const distance=radius/Math.sin(THREE.MathUtils.degToRad(camera.fov/2))/Math.min(camera.aspect,1)*1.04;
const dir={iso:new THREE.Vector3(-1.15,.95,1.25),top:new THREE.Vector3(0,1,.0001),front:new THREE.Vector3(-1,.0001,0),side:new THREE.Vector3(0,.0001,1)}[view];camera.position.copy(center).add(dir.normalize().multiplyScalar(distance));controls.target.copy(center);camera.near=Math.max(.001,radius/10000);camera.far=Math.max(100,distance*10);camera.updateProjectionMatrix();controls.update();}
function resize(){const w=stage.clientWidth,h=stage.clientHeight;renderer.setSize(w,h);camera.aspect=w/h;camera.updateProjectionMatrix();}
new ResizeObserver(resize).observe(stage);resize();
for(const b of document.querySelectorAll('[data-view]'))b.onclick=()=>fit(b.dataset.view);
document.querySelector('#fit').onclick=()=>fit();document.querySelector('#all').onclick=()=>{all();fit();};
document.querySelector('#source').onclick=()=>{for(const g of groups)g.object.visible=!g.id.startsWith('inferred_')&&!g.id.startsWith('restored_');updateCount();};
document.querySelector('#wire').onclick=e=>{wire=!wire;for(const m of meshes)m.material.wireframe=wire;e.currentTarget.setAttribute('aria-pressed',String(wire));};
document.querySelector('#panel-toggle').onclick=e=>{const open=document.body.classList.toggle('panel-open');e.currentTarget.setAttribute('aria-expanded',String(open));};
const ray=new THREE.Raycaster(),mouse=new THREE.Vector2();let down;
renderer.domElement.addEventListener('pointerdown',e=>down=[e.clientX,e.clientY]);
renderer.domElement.addEventListener('pointerup',e=>{if(!down||Math.hypot(e.clientX-down[0],e.clientY-down[1])>5)return;const rect=renderer.domElement.getBoundingClientRect();mouse.set((e.clientX-rect.left)/rect.width*2-1,-(e.clientY-rect.top)/rect.height*2+1);ray.setFromCamera(mouse,camera);
const shown=meshes.filter(m=>{let p=m;while(p){if(!p.visible)return false;p=p.parent;}return true;});const hit=ray.intersectObjects(shown,false)[0];clearSelection();if(hit){selected=hit.object;selected.material.emissive.setHex(0x37584d);document.querySelector('#selection strong').textContent=`${names[selected.userData.subsystem]} / ${selected.name||selected.parent.name||'Body'}`;}});
try{
meta=await fetch('assets/model.json').then(r=>{if(!r.ok)throw Error(`Metadata HTTP ${r.status}`);return r.json();});
const gltf=await new GLTFLoader().loadAsync('assets/spearhead-v2.glb',e=>{if(e.total)document.querySelector('#loading span').textContent=`Loading geometry · ${Math.round(e.loaded/e.total*100)}%`;});root=gltf.scene;
// Fusion glTF emits a primitive per face. Merge faces per named solid so
// the viewer renders ~316 body groups rather than issuing ~12,000 draw calls.
const mergeTargets=[];root.traverse(n=>{if(n.isGroup&&gltf.parser.associations.get(n)?.meshes!==undefined&&n.children.every(c=>c.isMesh))mergeTargets.push(n);});
for(const node of mergeTargets){const geometries=node.children.map(m=>{m.updateMatrix();const geo=m.geometry.clone().applyMatrix4(m.matrix);for(const name of Object.keys(geo.attributes))if(!['position','normal'].includes(name))geo.deleteAttribute(name);return geo;});const geometry=mergeGeometries(geometries,false);if(geometry){const mesh=new THREE.Mesh(geometry,node.children[0].material);mesh.name=node.name;node.clear();node.add(mesh);}for(const geo of geometries)geo.dispose();}
scene.add(root);
const groupHost=document.querySelector('#groups');groupHost.textContent='';
for(const entry of meta.groups){const object=root.getObjectByName(entry.id);if(!object)throw Error(`Missing subsystem: ${entry.id}`);const row=document.createElement('div');row.className='part-row';const label=document.createElement('label');const checkbox=document.createElement('input');checkbox.type='checkbox';checkbox.checked=true;checkbox.dataset.group=entry.id;label.append(checkbox,document.createTextNode(names[entry.id]));const qty=document.createElement('span');qty.className='qty';qty.textContent=entry.solids;const solo=document.createElement('button');solo.className='solo';solo.textContent='Solo';solo.title=`Isolate ${names[entry.id]}`;solo.setAttribute('aria-label',`Isolate ${names[entry.id]}`);row.append(label,qty,solo);groupHost.append(row);const group={...entry,object,checkbox};groups.push(group);checkbox.onchange=()=>{object.visible=checkbox.checked;updateCount();};solo.onclick=()=>{for(const g of groups)g.object.visible=g===group;updateCount();fit();};
object.traverse(node=>{if(node.isMesh){node.material=node.material.clone();node.material.color.set(colors[entry.id]);node.material.metalness=.05;node.material.roughness=.72;node.material.side=THREE.DoubleSide;node.userData.subsystem=entry.id;meshes.push(node);}});
}
fit();document.querySelector('#loading').hidden=true;document.querySelector('#loading').style.display='none';document.body.dataset.ready='true';
// Read-only diagnostics for repeatable delivery checks.
window.spearheadReview={get groups(){return groups.map(g=>({id:g.id,visible:g.object.visible,solids:g.solids}));},get meshes(){return meshes.length;},get bounds(){return new THREE.Box3().setFromObject(root).getSize(new THREE.Vector3()).toArray();},get camera(){return camera.position.toArray();},get target(){return controls.target.toArray();},get wire(){return wire;}};
}catch(error){document.querySelector('#loading').innerHTML='<strong>Could not load the model</strong><span>Reload the page to retry.</span>';console.error(error);}
function animate(){requestAnimationFrame(animate);controls.update();renderer.render(scene,camera);}animate();
