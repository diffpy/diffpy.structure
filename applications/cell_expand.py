"""This module includes functions to expand the cell structure to unitcell according to specified space group using cctbx package. Three functions are defined. These functions return a structure containing the  resulting structure and a string of constraints.
"""
import copy
import sys
import numpy as num
import numpy.linalg as numalg
from Structure import Structure
from Structure import Lattice
from Structure import Atom 
from cctbx import xray
from cctbx import crystal
from cctbx.array_family import flex

#############################################
    
def  super_structure_expand(spcgr, parameter_lists, parameter_initvalue_lists, init_structure=None, l=1, m=1, n=1):
    """ This function expand the cell using the space group given and   
    expand to a supercell using given l,m,n, it doesn't change the unitcell    
    for the lattice. Constraints are return in one string.     
    """
    if not isinstance(init_structure, Structure):
        raise RuntimeError, "expected instance of Structure"
    else:   
        ini_stru=init_structure
        ini_lat=init_structure.lattice
        unitcell=[float(y) for y in [ini_lat.a, ini_lat.b, ini_lat.c, ini_lat.alpha, ini_lat.beta, ini_lat.gamma] ]
        osymm = crystal.symmetry(unit_cell=unitcell, space_group_symbol=spcgr)
      
        scatterers=flex.xray_scatterer([xray.scatterer(label='H', u=0.2)]*len(ini_stru))

        constr_atom_para=['']*10
        constr_allatoms_para=[constr_atom_para]*len(ini_stru)
        atom_para_flag=[1]*10
        allatoms_para_flag=[atom_para_flag]*len(ini_stru)
        for atom_count in range(len(ini_stru)):
            for parameter_count in range(10):
                if   parameter_lists[atom_count][parameter_count].strip(' @')=='':
                    constr_allatoms_para[atom_count][parameter_count]=str(parameter_initvalue_lists[atom_count][parameter_count])
                    if parameter_initvalue_lists[atom_count][parameter_count]==0.0:
                        allatoms_para_flag[atom_count][parameter_count]=0                   
          
                else:
                    constr_allatoms_para[atom_count][parameter_count]=parameter_lists[atom_count][parameter_count].strip(' ')
      

      
        for atom in ini_stru: 
            print atom.name
            scatterers[ini_stru.index(atom)].label=atom.name
            scatterers[ini_stru.index(atom)].scattering_type=atom.element
            scatterers[ini_stru.index(atom)].site=(float(atom.xyz[0]),float(atom.xyz[1]),float(atom.xyz[2]))
            scatterers[ini_stru.index(atom)].u_iso=num.trace(atom.U)/3.0    
            scatterers[ini_stru.index(atom)].occupancy=float(atom.occupancy)
        
        ostru = xray.structure( special_position_settings=crystal.special_position_settings(crystal_symmetry=osymm), scatterers=scatterers )

# let's show the summary of the file anyway

        print "Summary of the current structure ..."
        ostru.show_summary().show_scatterers()
    

# step 2: write the scatterers (atoms)

#        j=10

        unitcell_atoms=[]   # store the information for atoms
        list_conts=[]
        origin_atom_count=-1
        unitcell_num_atoms=0
        for ascat in scatterers:
            origin_atom_count=origin_atom_count+1
            oequiv = ostru.sym_equiv_sites(ascat.site)
            numequiv = len(oequiv.coordinates())
            mtx_or=oequiv.original_site()
            for i in range(numequiv):
                New_atom=copy.deepcopy(ini_stru[origin_atom_count])
                one_op=oequiv.sym_op(i)
                mtx_ar=one_op.as_double_array()  # it is a 12 field tuple
      
      # the current approach is to calculate by ourself, though we can get those numbers from the oequiv too.
                x=mtx_or[0]*mtx_ar[0] + mtx_or[1]*mtx_ar[1] + mtx_or[2]*mtx_ar[2] + mtx_ar[9]
                y=mtx_or[0]*mtx_ar[3] + mtx_or[1]*mtx_ar[4] + mtx_or[2]*mtx_ar[5] + mtx_ar[10]
                z=mtx_or[0]*mtx_ar[6] + mtx_or[1]*mtx_ar[7] + mtx_or[2]*mtx_ar[8] + mtx_ar[11]
                New_atom.xyz=num.array([x,y,z], dtype=num.float64)
                unitcell_atoms.append(New_atom)
                unitcell_num_atoms = unitcell_num_atoms + 1
      #

                xstr=("".join([("%+f" % mtx_ar[ij]).rstrip("0") + '*(' + constr_allatoms_para[origin_atom_count][ij] +')' for ij in range(0,3) if (mtx_ar[ij] != 0.0 and allatoms_para_flag[origin_atom_count][ij] !=0) ])).lstrip('+')
                if mtx_ar[9] != 0.0 :
                    xstr=xstr + ("%+f" % mtx_ar[9]).rstrip('0')
                    
                ystr=("".join([("%+f" % mtx_ar[ij]).rstrip("0") + '*(' + constr_allatoms_para[origin_atom_count][ij - 3] +')' for ij in range(3,6) if (mtx_ar[ij] != 0.0 and allatoms_para_flag[origin_atom_count][ij - 3] !=0) ])).lstrip('+')
                if mtx_ar[10] != 0.0 :
                    ystr=ystr + ("%+f" % mtx_ar[10]).rstrip('0')

                zstr=("".join([("%+f" % mtx_ar[ij]).rstrip("0") + '*(' + constr_allatoms_para[origin_atom_count][ij - 6] +')' for ij in range(6,9) if (mtx_ar[ij] != 0.0 and allatoms_para_flag[origin_atom_count][ij - 6] !=0) ])).lstrip('+')
                if mtx_ar[11] != 0.0 :
                    zstr=zstr + ("%+f" % mtx_ar[11]).rstrip('0')

                ustr=constr_allatoms_para[origin_atom_count][3]
#               list_conts.extend([xstr+"\n", ystr+"\n", zstr+"\n", ustr + "\n"])
                list_conts.extend([xstr, ystr, zstr, ustr])


        allatoms=[] 
        allconstrain=[]
        sl=str(float(l))
        sm=str(float(m))
        sn=str(float(n))
        for ln in range(l):
            for mn in range(m):
                for nn in range(n):
                    move=num.array([ln, mn, nn], dtype=num.float64)
                    sln=''
                    if ln>0: sln=sln+'+'+str(float(ln)) 
                    smn=''
                    if mn>0: smn=smn+'+'+str(float(mn))
                    snn=''
                    if nn>0: snn=snn+'+'+str(float(nn))
       
                    for atom in unitcell_atoms:
                        moved_atom=copy.deepcopy(atom)
                        moved_atom.xyz=move+moved_atom.xyz
                        allatoms.append(moved_atom)
                        constrain_count=unitcell_atoms.index(atom)*4
                        allconstrain.extend([list_conts[constrain_count]+sln, list_conts[constrain_count+1]+smn, list_conts[constrain_count+2]+snn, list_conts[constrain_count+3]]) 
        super_structure=Structure(allatoms, ini_lat)
        super_structure.allconstrain=allconstrain
        return super_structure 
##########################################
  
def  cell_spcgr_expand(spcgr, parameter_lists, parameter_initvalue_lists, init_structure=None):
    """ This function only expand the structure using the space  
    group,constraints are returned. """ 
    if not isinstance(init_structure, Structure):
        raise RuntimeError, "expected instance of Structure"
    else:   
        ini_stru=init_structure
        ini_lat=init_structure.lattice
        unitcell=[float(y) for y in [ini_lat.a, ini_lat.b, ini_lat.c, ini_lat.alpha, ini_lat.beta, ini_lat.gamma] ]
        osymm = crystal.symmetry(unit_cell=unitcell, space_group_symbol=spcgr)
      
        scatterers=flex.xray_scatterer([xray.scatterer(label='H', u=0.2)]*len(ini_stru))
        constr_atom_para=['']*10
        constr_allatoms_para=[constr_atom_para]*len(ini_stru)
        atom_para_flag=[1]*10
        allatoms_para_flag=[atom_para_flag]*len(ini_stru)
        for atom_count in range(len(ini_stru)):
            for parameter_count in range(10):
                if  parameter_lists[atom_count][parameter_count].strip(' @')=='':
                    constr_allatoms_para[atom_count][parameter_count]=str(parameter_initvalue_lists[atom_count][parameter_count])
                    if  parameter_initvalue_lists[atom_count][parameter_count]==0.0:
                        allatoms_para_flag[atom_count][parameter_count]=0                   
          
                else:
                    constr_allatoms_para[atom_count][parameter_count]=parameter_lists[atom_count][parameter_count].strip(' ')
     
        for atom in ini_stru: 
            print atom.name
            scatterers[ini_stru.index(atom)].label=atom.name
            scatterers[ini_stru.index(atom)].scattering_type=atom.element
            scatterers[ini_stru.index(atom)].site=(float(atom.xyz[0]),float(atom.xyz[1]),float(atom.xyz[2]))
            scatterers[ini_stru.index(atom)].u_iso=num.trace(atom.U)/3.0    
            scatterers[ini_stru.index(atom)].occupancy=float(atom.occupancy)
        
        ostru = xray.structure( special_position_settings=crystal.special_position_settings(crystal_symmetry=osymm), scatterers=scatterers )

        print "Summary of the current structure ..."
        ostru.show_summary().show_scatterers()
    
# step 2: write the scatterers (atoms)


        unitcell_atoms=[]   # store the information for atoms
        list_conts=[]
        origin_atom_count=-1
        unitcell_num_atoms=0
        for ascat in scatterers:
            origin_atom_count=origin_atom_count+1
            oequiv = ostru.sym_equiv_sites(ascat.site)
            numequiv = len(oequiv.coordinates())
            mtx_or=oequiv.original_site()
            for i in range(numequiv):
                New_atom=copy.deepcopy(ini_stru[origin_atom_count])
                one_op=oequiv.sym_op(i)
                mtx_ar=one_op.as_double_array()  # it is a 12 field tuple
      
      # the currenet approach is to calculate by ourself, though we can get those numbers from the oequiv too.
                x=mtx_or[0]*mtx_ar[0] + mtx_or[1]*mtx_ar[1] + mtx_or[2]*mtx_ar[2] + mtx_ar[9]
                y=mtx_or[0]*mtx_ar[3] + mtx_or[1]*mtx_ar[4] + mtx_or[2]*mtx_ar[5] + mtx_ar[10]
                z=mtx_or[0]*mtx_ar[6] + mtx_or[1]*mtx_ar[7] + mtx_or[2]*mtx_ar[8] + mtx_ar[11]
                New_atom.xyz=num.array([x,y,z], dtype=num.float64)
                unitcell_atoms.append(New_atom)
                unitcell_num_atoms = unitcell_num_atoms + 1
      #
                xstr=("".join([("%+f" % mtx_ar[ij]).rstrip("0") + '*(' + constr_allatoms_para[origin_atom_count][ij] +')' for ij in range(0,3) if (mtx_ar[ij] != 0.0 and allatoms_para_flag[origin_atom_count][ij] !=0) ])).lstrip('+')
                if mtx_ar[9] != 0.0 :
                    xstr=xstr + ("%+f" % mtx_ar[9]).rstrip('0')
                    
                ystr=("".join([("%+f" % mtx_ar[ij]).rstrip("0") + '*(' + constr_allatoms_para[origin_atom_count][ij - 3] +')' for ij in range(3,6) if (mtx_ar[ij] != 0.0 and allatoms_para_flag[origin_atom_count][ij - 3] !=0) ])).lstrip('+')
                if mtx_ar[10] != 0.0 :
                    ystr=ystr + ("%+f" % mtx_ar[10]).rstrip('0')
                
                zstr=("".join([("%+f" % mtx_ar[ij]).rstrip("0") + '*(' + constr_allatoms_para[origin_atom_count][ij - 6] +')' for ij in range(6,9) if (mtx_ar[ij] != 0.0 and allatoms_para_flag[origin_atom_count][ij - 6] !=0) ])).lstrip('+')
                if mtx_ar[11] != 0.0 :
                    zstr=zstr + ("%+f" % mtx_ar[11]).rstrip('0')
                
                ustr=constr_allatoms_para[origin_atom_count][3]
                list_conts.extend([xstr+"\n", ystr+"\n", zstr+"\n", ustr + "\n"])
                list_conts.extend([xstr, ystr, zstr, ustr])


      # 
        unitcell=Structure(unitcell_atoms, ini_lat)
        unitcell.allconstrain=list_conts
        return unitcell     

##############################################   
def  super_unitcell_expand(unitcell_structure=None, l=1, m=1, n=1):
    """ This function create the new super-structure and thus changes the unitcell for lattice, constraints are returned in the list of strings. """
    if not isinstance(unitcell_structure, Structure):
        raise RuntimeError, "expected instance of Structure"
    else:   
        ini_stru=unitcell_structure
        ini_lat=unitcell_structure.lattice
        ini_const=unitcell_structure.allconstrain
        super_lattice=Lattice(ini_lat.a*l, ini_lat.b*m, ini_lat.c*n, ini_lat.alpha, ini_lat.beta, ini_lat.gamma)

    allatoms=[] 
    allconstrain=[]
    sl=str(float(l))
    sm=str(float(m))
    sn=str(float(n))
    for ln in range(l):
        for mn in range(m):
            for nn in range(n):
                fl=float(l)
                fm=float(m)
                fn=float(n)
                sfl=str(fl)
                sfm=str(fm)
                sfn=str(fn)
                sln=''
                if ln>0: sln=sln+'+'+str(float(ln)/fl) 
                smn=''
                if mn>0: smn=smn+'+'+str(float(mn)/fm)
                snn=''
                if nn>0: snn=snn+'+'+str(float(nn)/fn)  
                move=num.array([ln, mn, nn], dtype=num.float64)

                for atom in ini_stru:
                    moved_atom=copy.deepcopy(atom)
                    moved_atom.xyz=move+moved_atom.xyz
     #scale now
                    moved_atom.xyz[0]=moved_atom.xyz[0]/float(l)
                    moved_atom.xyz[1]=moved_atom.xyz[1]/float(m)
                    moved_atom.xyz[2]=moved_atom.xyz[2]/float(n)
           
                    allatoms.append(moved_atom)
                    constrain_count=ini_stru.index(atom)*4
                    allconstrain.extend(['('+ini_const[constrain_count]+')/'+sfl+sln, '('+ ini_const[constrain_count+1]+')/'+sfm+smn, '('+ini_const[constrain_count+2]+')/'+sfn+snn, ini_const[constrain_count+3]]) 

        super_unitcell=Structure(allatoms, super_lattice)
        super_unitcell.allconstrain=allconstrain  
        return super_unitcell      
#############################################
    
if __name__ == '__main__':
    print "OK"
