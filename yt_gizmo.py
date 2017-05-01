import yt
import pygadgetreader as pgt
import numpy as np
from yt.units import kiloparsec, Msun, kilometer, second, Gyr
from yt.units.yt_array import YTArray
from yt.utilities.math_utils import modify_reference_frame
from yt.utilities.cosmology import \
    Cosmology
from yt.data_objects.particle_filters import add_particle_filter

## change these as required

halo_particles = 553488
temp_unit = 98.54
numb_unit = 331.89365730278786 # hydrogen number density unit conv

#################

def read_snapshot(fileloc,field,ptype):
    data = pgt.readsnap(fileloc,field,ptype)
    print data, ptype, field, "size", len(data)
    return data


def gas_filter(pfilter,data):
	filter = np.logical_and(data["particle_type"] == 2, data["particle_type"] != None)
	return filter

def halo_gas_filter(pfilter,data):
	filter = np.logical_and(data["particle_type"] == 2, data["particle_ID"] <= halo_particles)
	return filter

def disk_gas_filter(pfilter,data):
	filter = np.logical_and(data["particle_type"] == 2, data["particle_ID"] > halo_particles)
	return filter

def star_filter(pfilter,data):
	filter = np.logical_or(data["particle_type"] == 1, data["particle_type"] == 4)
	return filter

def baryon_filter(pfilter,data):
	filter = np.logical_and(data["particle_type"] > 0, data["particle_type"] != None)
	return filter


def bulge_star_filter(pfilter,data):
	filter = np.logical_and(data["particle_type"] == 1, data["particle_type"] != None)
	return filter

def disk_star_filter(pfilter,data):
	filter = np.logical_and(data["particle_type"] == 4, data["particle_type"] != None)
	return filter

def dark_filter(pfilter,data):
	filter = np.logical_and(data["particle_type"] == 0, data["particle_type"] != None)
	return filter


def yt_gizmo(load_dir,snapno,halo_particles,read_list):
    # mass = Msun * 10^10
    # length = kpc
    # velocity = km/s
    # density = Msun * 10^10 / (kpc^3)

    alpha_thing = YTArray(6.5 * (10**34),"kg*m**(-3/2)*K**(-3/2)")
    # see alpha_2 from here http://astro.physics.uiowa.edu/~rlm/mathcad/addendum%204%20chap%2017%20stellar%20evolution%201.htm
    # see http://www.astronomy.ohio-state.edu/~ryden/ast825/ch1-3.pdf and http://astro.physics.uiowa.edu/~rlm/mathcad/addendum%204%20chap%2017%20stellar%20evolution%201.htm
    print "reading header"

    fileloc = load_dir + "/snapshot_" + snapno
    head = pgt.readheader(fileloc,"header")

    tot_part = head["npartThisFile"][0] + head["npartThisFile"][1] + head["npartThisFile"][2] + head["npartThisFile"][3]
    all_part = 0
    part_list = []
    part_read_list = []
    cum_part_list_right = []
    cum_part_list_left = []
    if "gas" in read_list:
        all_part += head["npartThisFile"][0]
        part_list.append(head["npartThisFile"][0])
        part_read_list.append("gas")
    if "dark" in read_list:
        all_part += head["npartThisFile"][1]
        part_list.append(head["npartThisFile"][1])
        part_read_list.append("dm")
    if "disk" in read_list:
        all_part += head["npartThisFile"][2]
        part_list.append(head["npartThisFile"][2])
        part_read_list.append("disk")
    if "bulge" in read_list:
        all_part += head["npartThisFile"][3]
        part_list.append(head["npartThisFile"][3])
        part_read_list.append("bulge")

    data_array = np.zeros([all_part,18])

    # make the part lists
    for j in range(0, len(part_list)):
        if j == 0:
            cum_part_list_left.append(0)
            cum_part_list_right.append(int(part_list[0]))
        if j == 1:
            cum_part_list_left.append(int(part_list[0]))
            cum_part_list_right.append(int(sum(part_list[:j+1])))
        if j > 1:
            cum_part_list_left.append(int(sum(part_list[:j])))
            cum_part_list_right.append(int(sum(part_list[:j+1])))

    print "reading snap"
    # TODO parallelize
    for i, ptype in enumerate(part_read_list):
        if ptype == "gas":
            # gas
            print "reading gas"
            data_array[cum_part_list_left[i]:cum_part_list_right[i],0] = 2    
        if ptype == "dm":
            # dark
            print "reading dark matter"
            data_array[cum_part_list_left[i]:cum_part_list_right[i],0] = 0  
        if ptype == "bulge":
            # bulge
            print "reading bulge stars"
            data_array[cum_part_list_left[i]:cum_part_list_right[i],0] = 1
        if ptype == "disk":
            # disk
            print "reading disk stars"
            data_array[cum_part_list_left[i]:cum_part_list_right[i],0] = 4

        # position

        data_array[cum_part_list_left[i]:cum_part_list_right[i],1:4] = read_snapshot(fileloc,"pos",ptype)
        # velocity
        data_array[cum_part_list_left[i]:cum_part_list_right[i],4:7] = read_snapshot(fileloc,"vel",ptype)
        # mass
        data_array[cum_part_list_left[i]:cum_part_list_right[i],7] = read_snapshot(fileloc,"mass",ptype) *  1e10
        # index
        data_array[cum_part_list_left[i]:cum_part_list_right[i],8] = read_snapshot(fileloc,"pid",ptype)

	data_array[cum_part_list_left[i]:cum_part_list_right[i],9] = read_snapshot(fileloc,"POT",ptype)

        if ptype == "gas":
            data_array[cum_part_list_left[i]:cum_part_list_right[i],10] = read_snapshot(fileloc,"u",ptype)
            data_array[cum_part_list_left[i]:cum_part_list_right[i],11] = read_snapshot(fileloc,"rho",ptype) *  1e10
            data_array[cum_part_list_left[i]:cum_part_list_right[i],12] = np.zeros(cum_part_list_right[i] - cum_part_list_left[i])
            data_array[cum_part_list_left[i]:cum_part_list_right[i],13] = read_snapshot(fileloc,"ne",ptype)
            data_array[cum_part_list_left[i]:cum_part_list_right[i],14] = read_snapshot(fileloc,"nh",ptype)
            data_array[cum_part_list_left[i]:cum_part_list_right[i],15] = data_array[cum_part_list_left[i]:cum_part_list_right[i],10] * temp_unit
            data_array[cum_part_list_left[i]:cum_part_list_right[i],16] = data_array[cum_part_list_left[i]:cum_part_list_right[i],11] * numb_unit / 1e10

        data_array[cum_part_list_left[i]:cum_part_list_right[i],17] = read_snapshot(fileloc,"TSTP",ptype)

    # convert into YT dataset
    mass = YTArray(data_array[:,7],"Msun").convert_to_units("g")

    print "making YT dataset"
    sim_data_verb = {"particle_type": YTArray(data_array[:,0],"dimensionless"),
                "particle_position_x": YTArray(data_array[:,1],"kpc").convert_to_units("cm"),
                "particle_position_y": YTArray(data_array[:,2],"kpc").convert_to_units("cm"),
                "particle_position_z": YTArray(data_array[:,3],"kpc").convert_to_units("cm"),
                "particle_velocity_x": YTArray(data_array[:,4], "km/s").convert_to_units("cm/s"),
                "particle_velocity_y": YTArray(data_array[:,5], "km/s").convert_to_units("cm/s"),
                "particle_velocity_z": YTArray(data_array[:,6], "km/s").convert_to_units("cm/s"),
                "particle_mass": mass,
                "cell_mass": mass,
                "particle_ID":  YTArray(data_array[:,8],"dimensionless"),
		"particle_potential": YTArray(data_array[:,9],"(km/s)**2").convert_to_units("(cm/s)**2"),
                "particle_energy": YTArray(data_array[:,10],"(km/s)**2").convert_to_units("(cm/s)**2"),
                "particle_density": YTArray(data_array[:,11],"Msun/(kpc**3)").convert_to_units("g/(cm**3)"),
                "particle_h": YTArray(data_array[:,12],"dimensionless"),
                "particle_ne": YTArray(data_array[:,13],"dimensionless"),
                "particle_nh": YTArray(data_array[:,14],"dimensionless"),
                "particle_temperature": YTArray(data_array[:,15],"K"),
                "particle_H_density": YTArray(data_array[:,16],"cm**-3").convert_to_units("cm**-3"),
		"particle_timestep": YTArray(data_array[:,17],"Gyr").convert_to_units("s")
            }

    sim_data_verb["particle_kinetic_energy"] = 0.5 * (np.power(sim_data_verb["particle_velocity_x"],2.0) + np.power(sim_data_verb["particle_velocity_y"],2.0) + np.power(sim_data_verb["particle_velocity_z"],2.0))

    # see http://www.astronomy.ohio-state.edu/~ryden/ast825/ch1-3.pdf and http://astro.physics.uiowa.edu/~rlm/mathcad/addendum%204%20chap%2017%20stellar%20evolution%201.htm
    jeans_mass =  alpha_thing.in_units("g*cm**(-3/2)*K**(-3/2)") * (sim_data_verb["particle_temperature"] ** (3.0/2.0) ) * ((sim_data_verb["particle_H_density"] / (1.0)) **(-1.0/2.0))

    sim_data_verb["particle_jeans_mass"] = YTArray(jeans_mass,"g")

    # final things
    min_pos = min(sim_data_verb["particle_position_x"].min().v, sim_data_verb["particle_position_y"].min().v, sim_data_verb["particle_position_z"].min().v) 
    max_pos = max(sim_data_verb["particle_position_x"].max().v, sim_data_verb["particle_position_y"].max().v, sim_data_verb["particle_position_z"].max().v)

    min_pos_check = min(sim_data_verb["particle_position_x"].min(), sim_data_verb["particle_position_y"].min(), sim_data_verb["particle_position_z"].min()) 
    max_pos_check = max(sim_data_verb["particle_position_x"].max(), sim_data_verb["particle_position_y"].max(), sim_data_verb["particle_position_z"].max())


    pos_val = max(abs(min_pos),abs(max_pos))

    bbox = 1.0 * np.array([[-pos_val,pos_val],
        [-pos_val,pos_val],
        [-pos_val,pos_val]])


    ytsnap = yt.load_particles(sim_data_verb, length_unit="cm", mass_unit="g", time_unit="s", bbox=bbox, n_ref=512)

    # filters
    print "adding particle filters"
    if "gas" in read_list:
        add_particle_filter("Gas", function=gas_filter, filtered_type="all", requires=["particle_type"])
        ytsnap.add_particle_filter("Gas")	

        add_particle_filter("HaloGas", function=halo_gas_filter, filtered_type="all", requires=["particle_type"])
        ytsnap.add_particle_filter("HaloGas")	

        add_particle_filter("DiskGas", function=disk_gas_filter, filtered_type="all", requires=["particle_type"])
        ytsnap.add_particle_filter("DiskGas")


    if "bulge" in read_list:
        add_particle_filter("BulgeStar", function=bulge_star_filter, filtered_type="all", requires=["particle_type"])
        ytsnap.add_particle_filter("BulgeStar")	

    if "disk" in read_list:
        add_particle_filter("DiskStar", function=disk_star_filter, filtered_type="all", requires=["particle_type"])
        ytsnap.add_particle_filter("DiskStar")

    if "bulge" in read_list and "gas" in read_list and "disk" in read_list:
        add_particle_filter("Baryon", function=baryon_filter, filtered_type="all", requires=["particle_type"])
        ytsnap.add_particle_filter("Baryon")	

    if "dark" in read_list:
        add_particle_filter("DarkMatter", function=dark_filter, filtered_type="all", requires=["particle_type"])
        ytsnap.add_particle_filter("DarkMatter")

    return ytsnap

if __name__ == "__main__":
    halo_particles = 553488
    dir_name = "./my_data"
    snap_name = "012"
    ytsnap = yt_gizmo(dir_name,"012",halo_particles,["gas","dark","bulge","disk"])

    ad = ytsnap.all_data()

    print "halo gas data"

    print len(ad["HaloGas","particle_ID"])
    print ad["HaloGas","particle_ID"].max(), ad["HaloGas","particle_ID"].min()
    print ad["HaloGas","particle_temperature"].max(), ad["HaloGas","particle_temperature"].min()
    print ad["HaloGas","particle_position_x"].max().in_units("kpc")

    print "disk data"
    print len(ad["DiskGas","particle_ID"])
    print ad["DiskGas","particle_ID"].max(), ad["DiskGas","particle_ID"].min()

    print ad["DiskGas","particle_position_x"].max().in_units("kpc")
    print ad["DiskGas","particle_temperature"].max(), ad["DiskGas","particle_temperature"].min()
    print ad["DiskGas","particle_density"].max(), ad["DiskGas","particle_density"].min()
    
