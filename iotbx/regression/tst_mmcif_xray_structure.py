"""
iotbx.pdb.mmcif.xray_structure_from_atom_site against the hierarchy route,
hierarchy_from_atom_site(block).extract_xray_structure(): the same scatterers
(label, scattering type, site, U, flags, occupancy, fp, fdp, multiplicity),
compared in label order since the direct route keeps the row order, over the
fixtures and model files of tst_mmcif_hierarchy_parity.py.
"""
import os
import libtbx.load_env
from libtbx.test_utils import Exception_expected
import iotbx.cif
from iotbx.cif.builders import CifBuilderError
from iotbx.pdb import mmcif
from iotbx.pdb import tst_mmcif_hierarchy as fixtures
from iotbx.regression import tst_mmcif_hierarchy_parity as parity


def scatterer_rows(xray_structure):
  result = []
  for sc in xray_structure.scatterers():
    result.append((sc.label, sc.scattering_type, tuple(sc.site), sc.u_iso,
      tuple(sc.u_star), sc.flags.use_u_iso(), sc.flags.use_u_aniso(),
      sc.occupancy, sc.fp, sc.fdp, sc.multiplicity()))
  return sorted(result)


def compare_block(block, label):
  hierarchy = mmcif.hierarchy_from_atom_site(block)
  # a block without a cell gets the box extract_xray_structure would make
  crystal_symmetry = hierarchy.generate_crystal_symmetry(
    mmcif.crystal_symmetry_builder(block).crystal_symmetry)
  ref = hierarchy.extract_xray_structure(crystal_symmetry=crystal_symmetry)
  new = mmcif.xray_structure_from_atom_site(block,
    crystal_symmetry=crystal_symmetry)
  assert new.crystal_symmetry().is_similar_symmetry(ref.crystal_symmetry())
  if scatterer_rows(new) != scatterer_rows(ref):
    raise AssertionError("%s: scatterers differ" % label)
  return new


def exercise_fixtures():
  for name, text in parity.parity_fixtures():
    compare_block(fixtures.first_block(text), "fixture " + name)


def exercise_model_files():
  n_aniso = 0
  for module, rel in parity.PARITY_FILES:
    fn = os.path.join(libtbx.env.dist_path(module), rel)
    assert os.path.isfile(fn), fn
    model = iotbx.cif.reader(file_path=fn).model()
    new = compare_block(list(model.values())[0], os.path.basename(fn))
    n_aniso += new.use_u_aniso().count(True)
  assert n_aniso > 0


def exercise_errors():
  block = fixtures.first_block(
    fixtures.atom_site_header + fixtures.input_two_atoms)
  try:
    mmcif.xray_structure_from_atom_site(block)
  except CifBuilderError as e:
    assert str(e) == "crystal symmetry is missing"
  else:
    raise Exception_expected


def run():
  exercise_fixtures()
  exercise_model_files()
  exercise_errors()


if __name__ == "__main__":
  run()
  print("OK")
