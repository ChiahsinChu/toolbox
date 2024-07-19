from typing import Union, List
import re
import glob
import sys

from pybtex.database import parse_file, BibliographyData


def extract_citation_keys(fnames: List[str]):
    """
    Grep all citation keys used in the tex file

    fname: str
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
):
    if tex_files is None:
        tex_files = glob.glob("./*.tex")
    if isinstance(tex_files, str):
        tex_files = [tex_files]

    citation_keys = extract_citation_keys(tex_files)
    bib_data = parse_file(bib_in_file)

    new_bib_data = BibliographyData()
    for key in citation_keys:
        if key in bib_data.entries:
            new_bib_data.entries[key] = bib_data.entries[key]

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
    export(bib_in_file, tex_files, bib_out_file)
