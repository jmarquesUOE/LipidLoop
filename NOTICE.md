# NOTICE — third-party spectral libraries distributed with LipidLoop

This file is kept at the root of the LipidLoop repository (github.com/jmarquesUOE/LipidLoop) and is
also included inside the spectral library archive released with it, because the MIT licence requires
its copyright notice to travel with every copy and CC BY 4.0 requires the attribution to do the same.

LipidLoop's **code** is MIT, Copyright (c) 2026 Jair G. Marques and Alex von Kriegsheim, University
of Edinburgh; see `LICENSE`. The **spectral libraries** are a separate matter. They come from three
sources and are redistributed under the terms below. `PROVENANCE_*.md` and `README.md` in this archive record, per library, where
it came from and exactly what was changed.

--------------------------------------------------------------------------------
1. LipiDex 1.1 — MIT licence
--------------------------------------------------------------------------------
Source:  https://github.com/coongroup/LipiDex  (src/msp_files/, src/libraries/)
Citation: Hutchins, P. D.; Russell, J. D.; Coon, J. J. LipiDex: An Integrated Software Package for
          High-Confidence Lipid Identification. Cell Syst. 2018, 6, 621-625.e5.
          https://doi.org/10.1016/j.cels.2018.03.011

Redistributed verbatim:
  LipiDex_HCD_Formic.msp, LipiDex_HCD_Acetate.msp, LipiDex_HCD_Hydroxy.msp,
  LipiDex_HCD_Plants.msp, LipiDex_HCD_ULCFA.msp, LipiDex_Splash_ISTD_Formic.msp,
  LipiDex_Splash_ISTD_Acetate.msp, LipidBlast_Formic.msp, LipidBlast_Acetate.msp
(Each is byte-identical to the file of the same name in coongroup/LipiDex at src/msp_files/.
LipidBlast_Formic.msp and LipidBlast_Acetate.msp are LipiDex's own in-silico libraries, generated
by LipiDex's LibraryGenerator from the templates in coongroup/LipiDex at src/libraries/; LipiDex 2
renames this library InSilico_LipidBlast. They are not the Fiehn Lab LipidBlast library file. The
in-silico approach is that of Kind, T. et al. LipidBlast in silico tandem mass spectrometry
database for lipid identification. Nat. Methods 2013, 10, 755-758.
https://doi.org/10.1038/nmeth.2551)

Modified copies derived from the above, distributed under the same licence:
  LipiDex_HCD_ULCFA_Deprotonated.msp   ([M-H]- subset of LipiDex_HCD_ULCFA.msp)
  Ceramides_Positive.msp               (positive-mode entries built from LipiDex's ceramide entries)
  DECOY_LipiDex_HCD_Acetate.msp, DECOY_LipidBlast_Acetate.msp,
  DECOY1_LipiDex_HCD_Formic.msp, DECOY1_LipiDex_HCD_Acetate.msp,
  DECOY1_LipiDex_HCD_Hydroxy.msp, DECOY1_LipiDex_HCD_ULCFA.msp,
  DECOY1_LipiDex_HCD_ULCFA_Deprotonated.msp, DECOY1_LipidBlast_Formic.msp,
  DECOY1_LipidBlast_Acetate.msp, DECOY2_HeadGroup_LipiDex_HCD_Formic.msp,
  DECOY2_HeadGroup_LipiDex_HCD_Acetate.msp

MIT License

Copyright (c) 2017 phutch89

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

--------------------------------------------------------------------------------
2. MS-DIAL Tandem Mass Spectral Atlas VS69 — CC BY 4.0
--------------------------------------------------------------------------------
Source:  Tsugawa, Hiroshi. "Lipidblast files". Zenodo, 10 April 2024.
         https://doi.org/10.5281/zenodo.10953284
         Files: MSDIAL-TandemMassSpectralAtlas-VS69-Pos.msp and -Neg.msp
Licence: Creative Commons Attribution 4.0 International (CC BY 4.0)
         https://creativecommons.org/licenses/by/4.0/
Please cite: Tsugawa, H.; Ikeda, K.; Takahashi, M.; Satoh, A.; Mori, Y.; Uchino, H.; Okahashi, N.;
         Yamada, Y.; Tada, I.; Bonini, P.; et al. A lipidome atlas in MS-DIAL 4.
         Nat. Biotechnol. 2020, 38, 1159-1163. https://doi.org/10.1038/s41587-020-0531-2

The following files contain material extracted from that atlas and MODIFIED. They are distributed
under CC BY 4.0.

  BileAcids_Negative.msp, BileAcids_Negative_Acetate.msp
    Changes: filtered to COMPOUNDCLASS BileAcid and BASulfate (47 of 792,757 entries); the name
    prefix "ST " relabelled "BA "; 11 entries with fewer than 3 peaks dropped; split by adduct into
    [M-H]- (25) and [M+CH3COO]- (11).

  UltraLongCer_EOS_Positive.msp, UltraLongHexCer_EOS_Positive.msp, AcylSM_Positive.msp,
  AcylHexCer_Positive.msp
    Changes: filtered to COMPOUNDCLASS Cer_EOS, HexCer_EOS, ASM and AHexCer (28,154 of 1.06 M
    entries); the subclass written into the compound name (Cer[EOS], HexCer[EOS], ASM); oxygen
    count rewritten ";2O" to ";O2"; the chain separator "/" rewritten "_"; the parenthesised
    esterified chain, "(O-14:0)" leading and "(FA 14:0)" trailing, rewritten as a chain.

  DECOY1_UltraLongCer_EOS_Positive.msp, DECOY1_UltraLongHexCer_EOS_Positive.msp,
  DECOY1_AcylSM_Positive.msp, DECOY1_AcylHexCer_Positive.msp
    Changes: target-decoy versions of the four libraries above, in which every chain-specific
    fragment mass is displaced by a whole number of methylene units while the precursor mass is
    kept. They are deliberately incorrect spectra, for false-discovery estimation only, and must
    not be used as reference spectra.

The full statement of changes is scripts/import_msdial_class.py and scripts/make_decoy_library.py
in the LipidLoop repository, together with PROVENANCE_BileAcids.md and PROVENANCE_UltraLong_Acyl.md
in this archive. No endorsement by the licensor is implied. The material is provided as-is; see the
CC BY 4.0 disclaimer of warranties at the URL above.

--------------------------------------------------------------------------------
3. Generated for LipidLoop — MIT licence
--------------------------------------------------------------------------------
Copyright (c) 2026 Jair G. Marques and Alex von Kriegsheim, University of Edinburgh, under the MIT
licence in the LipidLoop repository.

  Acylcarnitines_Positive.msp, Diacylglycerols_Positive.msp, NAcylAmides_Positive.msp,
  NAcylAmides_Negative.msp, OxTG_Positive.msp, Oxylipins_Negative.msp, Sterols_InSilico.msp,
  FreeFattyAcids_Negative.msp, Ganglioside_Negative_SumComposition.msp,
  and their DECOY1_ counterparts.

Ganglioside_Negative_SumComposition.msp is built by scripts/build_ganglioside_library.py from the
published Svennerholm glycan structures and elemental formulae, enumerating sphingoid bases C14-C22
against N-acyl chains C10-C32, with the three product ions measured in this laboratory's own
spectra. It reads no third-party file. It replaces an earlier library of the same name that took
its species list from LipiDex 2; that earlier file is not distributed. The rebuild also corrects an
inherited mass error in GD2-NGNA, which was 84.02 Da light.

Generated by the scripts/build_*_library.py and scripts/make_decoy_library.py scripts in the
LipidLoop repository, which are the reproducible record. Several are in-silico or partly in-silico
libraries whose limitations are stated in the PROVENANCE_*.md files in this archive; read those
before relying on an identification from them.

--------------------------------------------------------------------------------
4. LipiDex 2 — no material from it is distributed
--------------------------------------------------------------------------------
Source:  https://github.com/coongroup/LipiDex-2

The LipiDex 2 repository carries no licence file, states no terms in its README, and its
distributed archive contains no licence or end-user agreement, so no right to redistribute its
libraries can be asserted.

**No file in this archive is taken from or derived from LipiDex 2.** The libraries that once were
(LipiDex2_Ganglioside.msp, LipiDex2_FAHFA.msp, LipiDex2_Ornithine.msp and the decoys built from
them) have been removed, and the ganglioside library was rebuilt from published structures as
described in section 3 above.

LipiDex 2 is cited in the accompanying manuscript as prior work, which needs no licence:
Anderson, B. J.; Brademan, D. R.; He, Y.; Overmyer, K. A.; Coon, J. J. LipiDex 2 Integrates MSn
Tree-Based Fragmentation Methods and Quality Control Modules to Improve Discovery Lipidomics.
Anal. Chem. 2024, 96, 6715-6723. https://doi.org/10.1021/acs.analchem.4c00359
