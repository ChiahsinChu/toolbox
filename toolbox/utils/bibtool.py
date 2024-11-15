# SPDX-License-Identifier: LGPL-3.0-or-later
import glob
import re
import sys
from typing import List, Union

from doi2bib.crossref import get_bib_from_doi
from pybtex.database import BibliographyData, parse_file, parse_string
from tqdm import tqdm


def extract_citation_keys(fnames: List[str]):
    """
    Grep all citation keys used in the tex file

    Parameters
    ----------
    fnames : str
        Path to the tex file
    """
    citation_keys = set()
    for fname in fnames:
        with open(fname, "r") as file:
            tex_content = file.read()
        _citation_keys = set(re.findall(r"\\cite{([^}]+)}", tex_content))
        _citation_keys = {
            key.strip() for keys in _citation_keys for key in keys.split(",")
        }
        citation_keys.update(_citation_keys)
    return citation_keys


def export(
    bib_in_file: str,
    tex_files: Union[List[str], str] = None,
    bib_out_file: str = "export-ref.bib",
    online: bool = True,
):
    """
    Export a subset of the bib file that is used in the tex file

    Parameters
    ----------
    bib_in_file : str
        Path to the original bib file
    tex_files : Union[List[str], str]
        Path to the tex file or list of tex files
    bib_out_file : str
        Path to the output bib file
    online : bool
        If True, it will try to get the bib entry from the DOI
    """
    if tex_files is None:
        tex_files = glob.glob("./*.tex")
    if isinstance(tex_files, str):
        tex_files = [tex_files]

    citation_keys = extract_citation_keys(tex_files)
    bib_data = parse_file(bib_in_file)

    new_bib_data = BibliographyData()
    for kw in tqdm(citation_keys):
        if kw in bib_data.entries:
            if online:
                try:
                    out = get_bib_from_doi(
                        bib_data.entries[kw].fields["doi"], abbrev_journal=True
                    )
                    obj = parse_string(out[1], bib_format="bibtex")
                    for tmp_kw in obj.entries:
                        print(tmp_kw)
                    new_bib_data.entries[kw] = obj.entries[tmp_kw]
                except:
                    new_bib_data.entries[kw] = bib_data.entries[kw]
                new_bib_data.entries[kw].fields["title"] = bib_data.entries[kw].fields[
                    "title"
                ]
                try:
                    new_bib_data.entries[kw].fields["journal"] = bib_data.entries[
                        kw
                    ].fields["journal"]
                except KeyError:
                    pass
            else:
                new_bib_data.entries[kw] = bib_data.entries[kw]

    with open(bib_out_file, "w") as new_bib_file:
        new_bib_file.write(new_bib_data.to_string("bibtex"))


if __name__ == "__main__":
    """
    bibtool -i <bib_in_file> -o <bib_out_file> -tex <tex_files>
    if no tex_files is provided, it will search for all tex files in the current directory
    """
    args = sys.argv[1:]
    bib_in_file = args[args.index("-i") + 1]
    bib_out_file = (
        args[args.index("-o") + 1] if "-o" in args else "export-" + bib_in_file
    )
    tex_files = args[args.index("-tex") + 1] if "-tex" in args else None
    online = "--online" in args
    export(bib_in_file, tex_files, bib_out_file, online)
