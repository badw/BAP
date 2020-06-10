from pymatgen import Structure
import numpy as np
from pymatgen.core.sites import PeriodicSite, Site
import pandas as pd
import argparse
from tabulate import tabulate

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-f', '--file', type=str, default='POSCAR',
                                        help='path to input file')
    parser.add_argument('-v', '--verbose', default=False, action='store_true',
                                        help='verbose output - no grouping (default=False)')
    parser.add_argument('-a', '--atoms', nargs='+', default=None, help='atoms to be considered (in order)')
    parser.add_argument('-p', '--point', nargs='+', default=None, help='point to be considered i.e. a vacancy (x,y,z)')
    parser.add_argument('-r', '--radius', type=int, default=4, 
                                        help='radius for the point to find atoms (default=4)')
    parser.add_argument('-x', '--excel', default=False, action='store_true',
                                        help='output as excel friendly format (default=False)')
    parser.add_argument('-d', '--decimal', type=int, default=2,
                                        help='decimal place for the angles (default=2)')

    args = parser.parse_args()
    struct = Structure.from_file(args.file)
    struct_sites = struct.as_dict()['sites']
    atoms = args.atoms
    
    def specific_atom_centric(structure,atom1,atom2,atom3,radius):
        point = args.point
        if not point == None:
            site_search = struct.get_sites_in_sphere(point,include_image=True,r=args.radius)
            sites_init = [i[0].coords for i in site_search if i[0].as_dict()['species'][0]['element'] == atom1]
            sites_init_frac = [i[0].frac_coords for i in site_search if i[0].as_dict()['species'][0]['element'] == atom1] 
        else:    
            sites_init = [np.array(x['xyz']) for x in struct.as_dict()['sites'] if x['species'][0]['element'] == atom1 ] 
            sites_init_frac = [np.array(x['abc']) for x in struct.as_dict()['sites'] if x['species'][0]['element'] == atom1 ] 
        total_data = []
        for i,site in enumerate(sites_init): 
            def second_site_search(site1,new_radius):
                second_sites = struct.get_sites_in_sphere(site1,include_image=True,r=new_radius)
                new_second_sites = []
                new_second_sites_frac = []
                for x in second_sites:
                    if x[0].as_dict()['species'][0]['element'] == atom2:
                        if not np.array_equal(x[0].coords,site1):
                            new_second_sites.append(x[0].coords)
                            new_second_sites_frac.append(x[0].frac_coords)
                            
                return(new_second_sites,new_second_sites_frac)
     
            new_second_sites,new_second_sites_frac = second_site_search(site,radius)
                
            if new_second_sites == []:
                for new_radius in np.linspace(float(radius+0.1),radius+5,100):
                    new_second_sites,new_second_sites_frac = second_site_search(site,new_radius)
                    if not new_second_sites == []:
                        break
    
            for j,second_site in enumerate(new_second_sites):
                def third_site_search(site1,site2,new_radius):
                    third_sites = struct.get_sites_in_sphere(site2,include_image=True,r=new_radius)
                    new_third_sites = []
                    new_third_sites_frac = []
                    for x in third_sites:
                        if x[0].as_dict()['species'][0]['element'] == atom3:
                            if not np.array_equal(x[0].coords,site1):
                                if not np.array_equal(x[0].coords,site2):
                                    new_third_sites.append(x[0].coords)
                                    new_third_sites_frac.append(x[0].frac_coords)     
                    return(new_third_sites,new_third_sites_frac)
                
                new_third_sites,new_third_sites_frac = third_site_search(site,second_site,radius)
                
                if new_third_sites == []:
                        for new_radius in np.linspace(float(radius+0.1),radius+5,100):
                            new_third_sites,new_third_sites_frac = third_site_search(site,second_site,new_radius)
                            if not new_third_sites == []:
                                break
    
                for k,third_site in enumerate(new_third_sites):
                    ba = site - second_site
                    bc = third_site - second_site
                    
                 
                    cosine_angle = np.dot(ba, bc) / (np.linalg.norm(ba) * np.linalg.norm(bc))
                    angle = np.around(np.degrees(np.arccos(cosine_angle)),decimals=args.decimal)
                    total_data.append([sites_init_frac[i],new_second_sites_frac[j],new_third_sites_frac[k],angle])
                        
    
                        
                    
                np.set_printoptions(formatter={'float': lambda x: "{0:0.3f}".format(x)}) 
      
    
        df = pd.DataFrame(total_data, columns = [atom1, atom2, atom3,'angle'])
        df = df.sort_values(by=['angle'])
        df = df.reset_index(drop=True)
        if args.verbose == False:
            df = df.drop_duplicates(subset='angle', keep='first')
            df = df.reset_index(drop=True)
            return(df)
        else:
            return(df)
    
    bond_angle_data = specific_atom_centric(struct,atoms[0],atoms[1],atoms[2],3)
    if args.excel == True:
        print(bond_angle_data.to_csv(index=True,header=True))
    else:
        print(tabulate(bond_angle_data, headers='keys', tablefmt='psql'))

if __name__ == "__main__":
    main()
