# %%
from openff.toolkit.typing.engines.smirnoff import forcefield
from pprint import pprint
pprint(forcefield.get_available_force_fields())
# %%
# 処理したいPDBファイルの指定
pdb_file = 'Data/3poz.pdb'

# 必要なライブラリのインポート
import pdbfixer

# PDBFixerでPDBファイルを読み込む
fixer = pdbfixer.PDBFixer(pdb_file)

# 不要な構造の削除
fixer.removeHeterogens()

# 欠けている残基のチェック（欠損原子の確認のためにも必要）
fixer.findMissingResidues()

# タンパク質末端の欠けている残基を取り除く処理
chains = list(fixer.topology.chains())
keys = fixer.missingResidues.keys()
for key in list(keys):
    chain = chains[key[0]]
    if key[1] == 0 or key[1] == len(list(chain.residues())):
        del fixer.missingResidues[key]

# 非標準な残基が含まれているか確認、あれば標準的なものに置き換える
fixer.findNonstandardResidues()
fixer.replaceNonstandardResidues()

# 欠けている原子の確認、あれば追加する
fixer.findMissingAtoms()
fixer.addMissingAtoms()

# 水素原子の付与（pHを設定する）
ph = 7.0
fixer.addMissingHydrogens(ph)

print(type(fixer))
# %%
from rdkit import Chem
from rdkit.Chem import AllChem

# PDBファイルを読み込んで残基ごとにスプリット
rdkit_pdb = Chem.MolFromPDBFile('Data/3poz.pdb')
rdkit_pdb_split = Chem.SplitMolByPDBResidues(rdkit_pdb)

# 残基 03Pの抜きだし
pdb_tak285 = rdkit_pdb_split["03P"]

Chem.Draw.MolToImage(pdb_tak285)
# %%
# 中途半端に残ってる水素原子を消しておく
pdb_tak285_noH = Chem.RemoveHs(pdb_tak285)

# SMILESからリファレンス構造を作成
smi_tak285 = "CC(C)(CC(=O)NCCn1ccc2c1c(ncn2)Nc3ccc(c(c3)Cl)Oc4cccc(c4)C(F)(F)F)O"
ref_tak285 = Chem.MolFromSmiles(smi_tak285)

# リファレンスをもとに結合次数をアサイン
prep_tak285 = AllChem.AssignBondOrdersFromTemplate(ref_tak285, pdb_tak285)
prep_tak285.AddConformer(pdb_tak285.GetConformer(0))

# 水素原子の付加
prep_tak285 = Chem.AddHs(prep_tak285, addCoords=True)

Chem.Draw.MolToImage(prep_tak285)
# %%
from openmm import *
from openmm.app import *

FF = ForceField('amber14-all.xml', 'amber14/tip3pfb.xml')
# %%
from openff.toolkit.topology import Molecule

off_tak285 = Molecule.from_rdkit(prep_tak285, allow_undefined_stereo=True, hydrogens_are_explicit=True)

print(type(off_tak285))
# %%
from openmmforcefields.generators import GAFFTemplateGenerator

gaff_temp = GAFFTemplateGenerator(off_tak285)
FF.registerTemplateGenerator(gaff_temp.generator)
# %%
tak285_template = gaff_temp.ge
# %%
print("version info : ", gaff_temp.gaff_version)
print("installed versions : ", GAFFTemplateGenerator.INSTALLED_FORCEFIELDS)
# %%
# 分子の名前（"LIG"）を追加
off_tak285.name = "LIG"

# 原子の名前を追加
# 辞書を用意しておいて、元素記号＋数字形式にする
element_counter_dict = {}
for off_atom, rdkit_atom in zip(off_tak285.atoms, prep_tak285.GetAtoms()):
    element = rdkit_atom.GetSymbol()
    if element in element_counter_dict.keys():
        element_counter_dict[element] += 1
    else:
        element_counter_dict[element] = 1
    off_atom.name = element + str(element_counter_dict[element])
# %%
# OpenMMの単位系導入
from openmm import unit

# トポロジー（topology）の取り出しと変換
off_top = off_tak285.to_topology()
omm_top = off_top.to_openmm()

# 座標（positions）の取り出し
omm_pos = off_tak285.conformers[0]

# 単位をÅからnmへ
for atom in omm_pos:
    coords = atom / atom.units
    atom = (coords / 10.0) * unit.nanometers

# TopologyとPositionsをまとめてOpenMMのオブジェクトに変換
omm_tak285 = app.Modeller(omm_top, omm_pos)
# %%
import mdtraj as md

# タンパク質とリガンドのトポロジーをそれぞれMdTrajに変換
md_protein_top = md.Topology.from_openmm(fixer.topology) 
md_ligand_top = md.Topology.from_openmm(omm_tak285.topology) 

# joinでひとつのトポロジーにまとめる
md_complex_top = md_protein_top.join(md_ligand_top)

# OpenMMオブジェクトにする
omm_complex_top = md_complex_top.to_openmm()
# %%
import numpy as np

# タンパクとリガンドを合わせた総数
total_atoms = len(fixer.positions) + len(omm_tak285.positions)

# positionsを格納するQuantity配列の準備
complex_positions = unit.Quantity(np.zeros([total_atoms, 3]), unit=unit.nanometers)

# タンパク質のpositionsを追加
complex_positions[: len(fixer.positions)] = fixer.positions
# 後に続けてリガンドのpositionsを追加
complex_positions[len(fixer.positions) :] = omm_tak285.positions
# %%
comp_model = app.Modeller(omm_complex_top, complex_positions)
# %%
# comp_model.addSolvent(FF, padding=1.0 * unit.nanometers, ionicStrength=0.15 * unit.molar)

# %%
# 処理後の状態（トポロジー、原子の位置）をPDBファイルで出力
top = comp_model.getTopology()
pos = comp_model.getPositions()
app.PDBFile.writeFile(top, pos, open('Data/gaff_complex_processed.pdb', 'w'))
# %%
from openmm import *
from openmm.app import *
from openmm.unit import *

# シミュレーションの準備
## ForceFiledからSystemを構築する
FF.registerTemplateGenerator(gaff_temp.generator)
system = FF.createSystem(comp_model.topology, nonbondedMethod=PME, nonbondedCutoff=1.0*nanometers, constraints = HBonds, rigidWater=True, ewaldErrorTolerance=0.0005)

## 圧力制御の機能をaddForceする
system.addForce(MonteCarloBarostat(1.0*atmospheres, 300*kelvin, 25))

## 積分計算の設定
integrator = LangevinMiddleIntegrator(300*kelvin, 1.0/picosecond, 0.002*picoseconds)
integrator.setConstraintTolerance(0.000001)

## 計算環境をPlatformで設定
platform = Platform.getPlatformByName('CUDA')
platformProperties = {'Precision': 'single'}

## simulationの構築
simulation = Simulation(comp_model.topology, system, integrator, platform, platformProperties)
simulation.context.setPositions(comp_model.positions)

# XML結果出力ための設定
with open("system.xml", mode="w") as file:
    file.write(XmlSerializer.serialize(system))
with open("integrator.xml", mode="w") as file:
    file.write(XmlSerializer.serialize(integrator))

# エネルギー最小化の実行
print('Performing energy minimization...')
simulation.minimizeEnergy()

# 平衡化のステップ
print('Equilibrating...')
simulation.context.setVelocitiesToTemperature(300*kelvin)
simulation.step(50000)

# シミュレーション本番の実行の前にreporterオブジェクトを追加する
print('Simulating...')

dcdReporter = DCDReporter('trajectory.dcd', 10000)
dataReporter = StateDataReporter('log.txt', 1000, totalSteps=500000, step=True, speed=True, progress=True, potentialEnergy=True, temperature=True, separator='\t')
checkpointReporter = CheckpointReporter('checkpoint.chk', 10000)

simulation.reporters.append(dcdReporter)
simulation.reporters.append(dataReporter)
simulation.reporters.append(checkpointReporter)
simulation.currentStep = 0

# シミュレーション本番の実行
simulation.step(500000)

# シミュレーション最後の状態をXMLで書き出す
simulation.saveState("final_state.xml")

# 最後の状態をPDBx/mmcifで出力する設定
state = simulation.context.getState(getPositions=True, enforcePeriodicBox=system.usesPeriodicBoundaryConditions())
with open("final_state.pdbx", mode="w") as file:
    PDBxFile.writeFile(simulation.topology, state.getPositions(), file)
# %%
from rdkit import Chem
from openff.toolkit.topology import Molecule

# Load your RDKit molecule (e.g., from a SMILES string)
rdkit_molecule = Chem.MolFromSmiles('CCO')  # Example: ethanol

# Convert RDKit molecule to OpenFF Molecule
off_molecule = Molecule.from_rdkit(rdkit_molecule)
from openff.toolkit.typing.engines.smirnoff import ForceField

# Load the OpenFF "Parsley" force field
forcefield = ForceField('openff_unconstrained-1.3.0.offxml')

# Create an OpenFF Topology
off_topology = off_molecule.to_topology()

# Create OpenMM system
omm_system = forcefield.create_openmm_system(off_topology)

# %%
from openmm import XmlSerializer

with open('molecule_system.xml', 'w') as f:
    f.write(XmlSerializer.serialize(omm_system))

# %%
from rdkit import Chem
from openff.toolkit.topology import Molecule

# Load your RDKit molecule (e.g., from a SMILES string)
rdkit_molecule = Chem.MolFromSmiles('CCO')  # Example: ethanol

# Generate 3D coordinates
Chem.AddHs(rdkit_molecule)
Chem.EmbedMolecule(rdkit_molecule)

# Convert RDKit molecule to OpenFF Molecule
off_molecule = Molecule.from_rdkit(rdkit_molecule)

# %%
from rdkit import Chem
from openff.toolkit.topology import Molecule

# Load your RDKit molecule (e.g., from a SMILES string)
rdkit_molecule = Chem.MolFromSmiles('CCO')  # Example: ethanol

# Generate 3D coordinates
Chem.AddHs(rdkit_molecule)
# Chem.EmbedMolecule(rdkit_molecule)

# Convert RDKit molecule to OpenFF Molecule
off_molecule = Molecule.from_rdkit(rdkit_molecule)

# %%
from openmmforcefields.generators import GAFFTemplateGenerator

gaff_template = GAFFTemplateGenerator(off_molecule)

# %%
from openmm.app import ForceField

forcefield = ForceField('amber14/protein.ff14SB.xml', 'amber14/tip3p.xml')
forcefield.registerTemplateGenerator(gaff_template.generator)

# %%
from openmm.app import PDBFile
from openmm import XmlSerializer

# Create a PDB file for the molecule
with open('molecule.pdb', 'w') as pdb_file:
    PDBFile.writeFile(off_molecule.to_topology().to_openmm(), off_molecule.conformers[0], pdb_file)

# Load the PDB file
pdb = PDBFile('molecule.pdb')

# Create the system
system = forcefield.createSystem(pdb.topology)

# Serialize the system to an XML file
with open('molecule_system.xml', 'w') as xml_file:
    xml_file.write(XmlSerializer.serialize(system))

# %%
