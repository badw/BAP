from pymatgen import Structure
import numpy as np
import argparse
from tabulate import tabulate
from pymatgen.core.sites import PeriodicSite

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-f', '--file', type=str, default='POSCAR',
                                        help='path to input file')
    parser.add_argument('-p', '--planar_atoms', nargs='+', default=None, help='atom indexes to be considered as the plane QRS (choose 3) (in order)')
    parser.add_argument('-a', '--projected_atom', type=str, default=None, help='atom index to be projected onto the plane')
    parser.add_argument('-d', '--decimal', type=int, default=2,
                                        help='decimal place for the angles (default=2)')
    parser.add_argument('-o', '--output', action='store_true', default=False,
                                        help='generate a new POSCAR file with the projected atom (as H)')

    args = parser.parse_args()
    struct = Structure.from_file(args.file)
    plane = [int(x) for x in args.planar_atoms]
    atoms = struct.sites[int(args.projected_atom)].frac_coords
    
    planaratoms = []
    species = []
    for i,at in enumerate(struct.sites):
        if i in plane:
            planaratoms.append(at.frac_coords)
            species.append(at.specie) 
    
    species.append(struct.sites[int(args.projected_atom)].specie)
    
    def angles(Q,R,S,O):
        RQ_vec = R-Q
        RS_vec = R-S
        RO_vec = R-O
        RQ_mag = np.linalg.norm(RQ_vec)
        RS_mag = np.linalg.norm(RS_vec)
        RO_mag = np.linalg.norm(RO_vec)

        n = np.cross(RQ_vec,RS_vec) / (RQ_mag * RS_mag)

        PO_mag = np.dot(np.array(n),RO_vec)
        PO_vec = PO_mag * np.array(n)
        RP_vec = RO_vec - PO_vec
        RP_mag = np.linalg.norm(RP_vec)

        P = R - RP_vec

        thetaz = np.dot(RQ_vec,RP_vec)/(RQ_mag * RP_mag)
        angle = np.arccos(thetaz)

        return(np.degrees(angle),P)

    tilt,projection =  angles(planaratoms[0],planaratoms[1],planaratoms[2],atoms) 
    
    table = [['Q','{}{}'.format(species[0],plane[0]),planaratoms[0]],
             ['R','{}{}'.format(species[1],plane[1]),planaratoms[1]],
             ['S','{}{}'.format(species[2],plane[2]),planaratoms[2]],
             ['original atom (O)','{}{}'.format(species[3],int(args.projected_atom)),atoms],
             ['projected atom (P)','{}{}'.format('-','-'),projection],
             ['angle','<(RQ,RP)',tilt]]
    
    print('''
    Q               
    |    P           
    |---/   O         
    |  /  >           
    | / >             
    |/>               
    R-----------S''')

    print(tabulate(table))
    
     
    if args.output == True:
        struct.append("H",projection)
        struct.to(filename='POSCAR_with_projection.vasp',fmt='poscar')
        print('projection outputted to file (P = Hydrogen)')

if __name__ == "__main__":
    main()
