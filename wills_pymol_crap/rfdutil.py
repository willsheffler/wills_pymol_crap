from pymol import cmd
from wills_pymol_crap.sym_util import aligncx
import yaml

def make_gp_flags(sele, path='./sym_gp'):
    path = os.path.abspath(path)
    pdbfile = f'{path}.pdb'
    conffile = f'{path}.yaml'
    nfold = 2
    conf = dict(defaults='base sym viz _self_'.split(),
                hydra=dict(searchpath='pkg://rf_diffusion/config/inference'),
                viz=dict(rfold_iter_end=False, diffusion_step=False),
                inference=dict(input_pdb=pdbfile, contig_as_guidepost=True),
                sym=dict(
                    symid=f'C{nfold}',
                    sympair_protein_only=False,
                    move_unsym_with_asu=False,
                    start_radius=0.0,
                    contig_is_symmetric=True,
                    fit=False,
                ),
                contigmap=dict(has_termini=[True] * nfold),
                diffuser=dict(T=20))
    add_tip_config(conf, sele)
    cmd.save(pdbfile)
    print(yaml.dump(conf))
    with open(conffile, 'w') as out:
        yaml.dump(conf, out)

def add_tip_config(conf, sele):
    mdl = cmd.get_model(f'({sele})')
    ligmdl = cmd.get_model(f'(({sele}) and het)')
    ligresn, ligresi = zip(*list({(a.resn, a.resi) for a in mdl.atom if a.hetatm}))
    if len(ligresn) != 1 or len(ligresi) != 1:
        raise ValueError(f'must be exactly one ligand {ligresn} {ligresi}')
    ligresn, ligresi = ligresn[0], ligresi[0]
    reslist = list(sorted({(a.chain, int(a.resi), a.resn) for a in mdl.atom if not a.hetatm}))
    conf['contigmap']['contigs'] = []
    conf['contigmap']['contig_atoms'] = {}
    conf['inference']['ligand'] = ligresn
    for chain, resi, resn in reslist:
        print(chain, resi, resn)
        if resn not in resn_tip_atoms:
            print(f'no tip atoms for resn {resn}')
            continue
        tip = resn_tip_atoms.get(resn, [])
        cmd.select('test',
                   f'(({sele}) and chain {chain} and resi {resi} and name {"+".join(tip)})',
                   merge=True)
        conf['contigmap']['contigs'].append(chain + str(resi))
        conf['contigmap']['contig_atoms'][chain + str(resi)] = ','.join(tip)
    conf['contigmap']['contigs'] = [','.join(conf['contigmap']['contigs'])]

resn_tip_atoms = dict(
    ARG='NE CZ NH1 NH2'.split(),
    ASN='CG OD1 ND2'.split(),
    GLN='CD OE1 NE2'.split(),
    ASP='CG OD1 OD2'.split(),
    GLU='CD OE1 OE2'.split(),
    LYS='CD CE NZ'.split(),
    SER='CB OG'.split(),
    # VAL=''.split(),
)
