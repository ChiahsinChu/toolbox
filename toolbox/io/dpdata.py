# SPDX-License-Identifier: LGPL-3.0-or-later
import glob
import os
from typing import List

import numpy as np
from tqdm import tqdm


def wannier_v2_to_v3(dname, sel_type: List[str]):
    fnames = glob.glob(os.path.join(dname, "**/atomic_dipole.npy"), recursive=True)
    fnames.sort()
    for fname in tqdm(fnames):
        dname = os.path.dirname(fname)
        # print(dname)
        _type_map = np.loadtxt(os.path.join(dname, "../type_map.raw"), dtype=str)
        atype = np.loadtxt(os.path.join(dname, "../type.raw"), dtype=int)
        symbols = _type_map[atype]
        sel_ids = []
        for t in sel_type:
            sel_ids.extend(np.where(symbols == t)[0].tolist())
        old_data = np.load(fname)
        n_frames = old_data.shape[0]
        new_data = np.zeros_like(np.load(os.path.join(dname, "coord.npy"))).reshape(
            [n_frames, -1, 3]
        )
        # np.copyto(new_data[:, sel_ids], old_data.reshape(n_frames, -1, 3))
        new_data[:, sel_ids] = old_data.reshape(n_frames, -1, 3)
        os.rename(
            os.path.join(dname, "atomic_dipole.npy"),
            os.path.join(dname, "v2_atomic_dipole.npy"),
        )
        np.save(
            os.path.join(dname, "atomic_dipole.npy"), new_data.reshape([n_frames, -1])
        )
