# SPDX-License-Identifier: LGPL-3.0-or-later
import unittest

from pathlib import Path

from toolbox.io.cp2k import Cp2kOutput


class TestCp2kInput(unittest.TestCase):
    pass


class Cp2kOutputTest:
    def test_energy(self):
        self.assertEqual(
            self.cp2k_out.energy,
            self.ref_data["energy"],
        )

    def test_multiplicity(self):
        self.assertEqual(
            self.cp2k_out.multiplicity,
            self.ref_data["multiplicity"],
        )

    def test_fermi(self):
        self.assertEqual(
            self.cp2k_out.fermi,
            self.ref_data["fermi"],
        )

    def test_uks(self):
        self.assertEqual(
            self.cp2k_out.uks,
            self.ref_data["uks"],
        )

    def test_scf_loop(self):
        self.assertEqual(
            self.cp2k_out.scf_loop,
            self.ref_data["scf_loop"],
        )
        
    def test_charge(self):
        self.assertEqual(
            self.cp2k_out.charge,
            self.ref_data["charge"],
        )


class TestCp2kCube(unittest.TestCase):
    pass


class TestCp2kHartreeCube(unittest.TestCase):
    pass


class TestRKSCp2kOutput(Cp2kOutputTest, unittest.TestCase):
    def setUp(self) -> None:
        fname = str(Path(__file__).parents[1] / "data" / "cp2k_v7.1_ot_rks" / "output")
        self.cp2k_out = Cp2kOutput(fname)
        self.ref_data = {
            "multiplicity": 1,
            "fermi": 1.326202,
            "uks": False,
            "energy": -6613.877097770478031 * 2.72113838565563e01,
            "scf_loop": 8,
            "charge": 0,
        }


class TestUKSCp2kOutput(Cp2kOutputTest, unittest.TestCase):
    def setUp(self) -> None:
        fname = str(
            Path(__file__).parents[1] / "data" / "cp2k_v2024.3_ot_uks" / "output"
        )
        self.cp2k_out = Cp2kOutput(fname)
        self.ref_data = {
            "multiplicity": 2,
            "fermi": 4.056627,
            "uks": True,
            "energy": -4172.293059921372333 * 2.72113838565563e01,
            "scf_loop": 1,
            "charge": 0,
        }

class TestRKSSmearCp2kOutput(Cp2kOutputTest, unittest.TestCase):
    def setUp(self) -> None:
        fname = str(
            Path(__file__).parents[1] / "data" / "cp2k_v2024.3_smear_rks" / "output"
        )
        self.cp2k_out = Cp2kOutput(fname)
        self.ref_data = {
            "multiplicity": 1,
            "fermi": 0.16599014083981 * 2.72113838565563e01,
            "uks": False,
            "energy": -4172.132809629093572 * 2.72113838565563e01,
            "scf_loop": 39,
            "charge": -1,
        }


if __name__ == "__main__":
    unittest.main()
