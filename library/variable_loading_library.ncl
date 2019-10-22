;
;  Universal library for loading variables from a list of netcdf files
;

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRFUserARW.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/contrib/time_axis_labels.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRF_contributed.ncl"

load "/Volumes/Odin/YRD_China/data_plotting/scripts/library/WRFCHEM_misc_stats.ncl"

begin

;;; useful physical constants
m_dry = 28.97           ; molar mass dry air (g/mol)
;m_no3 = 62.005          ; molar mass NO3
m_nh3 = 17.03           ; 
;m_nh4 = 18.04           ; etc...
;m_so4 = 96.06
;m_cl  = 35.453
m_vap = 18.01528     ; H2O
m_oxygen = 16.0
m_carbon = 12.0
m_h2o = m_vap
;m_na  = 22.9898

;;; use molar weights from MOSAIC, to ensure as much compatibility as possible for pH
m_no3 = 62.0
m_so4 = 96.0
m_cl  = 35.5
m_nh4 = 18.0
m_na  = 23.0


m_vbs = 250.0   ; use a single molar mass for VBS compounds - change this if WRF-Chem changes!!!


; Aerosol densities (from MOSAIC data file module_data_mosaic_asect). in g/cm^3
dens_so4_aer = 1.80
dens_no3_aer = 1.80
dens_cl_aer  = 2.20
; dens_msa_aer = 1.80
; dens_co3_aer = 2.60
dens_nh4_aer = 1.80
dens_na_aer  = 2.20
; dens_ca_aer  = 2.60
dens_oin_aer = 2.60
dens_oc_aer  = 1.00
dens_bc_aer  = 1.70
dens_water_aer  = 1.0


N_A = 6.02214e23    ; avogadro's number mol^-1
R = 8.3144621       ; gas constant J/molK


;; conversion factors

;;;;; kg-1 (dry air) to m-3 (sample volume - standard atmosphere)
; Using R_dry = 287.058 J kg^-1 K^-1
; T = 273.15K
; p = 101325 Pa = 101325 J m^-3

; V = M * (R_d*T)/p = M * 0.7738, V^-1 = 1.2923 * M^-1
; Therefore multiply by factor 1.2923 to convert "/kg_dry" to "/sm^3"  
conv_kg_sm3 = 1.2923



;-------------------------------------------------------------------
;     Convert from ppm to molecules cm^-3
; !!! WARNING - Must have pressure in Pa !!!

undef("ppm_moleconc_conv")
function ppm_moleconc_conv(p:numeric, Temp:numeric)	; pressure [Pa], temperature [K]
begin

	if (p@units .eq. "Pa")
		 ppm_conv_fac = 1e-12*N_A*p/(Temp*R)
	elseif (p@units .eq. "hPa")
		 ppm_conv_fac = 1e-10*N_A*p/(Temp*R)
	else
		print("ppm_conv: must have pressure in units Pa or hPa")
	end if		

	 return(ppm_conv_fac)

end


;-------------------------------------------------------------------;
; get_grid_thickness	-	return the vertical thickness(del_z) of grid cells (in m)
; SAN - 19/02/14 - Reframed to be calculated explicitly from geopotential at grid 
; boundaries (unstaggered)

; Height = full geopotential / 9.81

undef("get_grid_thickness")
function get_grid_thickness(PH:numeric, PHB:numeric, height:numeric)
local z_stag, del_z, dims, nk, iz
begin

	; calculate the height at the staggered grid boundaries
	z_stag = (PH + PHB)/9.81

	; get dimension sizes of the domain
	dims = dimsizes(height)

	; create array for model thickness same size/type as the unstaggered height array
	del_z = new(dims,typeof(height))
	
	
	; check if the variable has 3 or 4 dimensions (i.e. does it have a time dimension or not)
	if(dimsizes(dims).eq.4)
		nk = dims(1)

		; loop through the vertical, calculating thickness at each level:
		do iz = 0, nk-1
			del_z(:, iz, :, :) = z_stag(:, iz+1, :, :) - z_stag(:, iz, :, :)
		end do
	
	elseif (dimsizes(dims).eq.3)
		nk = dims(0)

		; loop through the vertical, calculating thickness at each level:
		do iz = 0, nk-1
			del_z(iz, :, :) = z_stag(iz+1, :, :) - z_stag(iz, :, :)
		end do

	end if
	

	del_z@units = "m"
	del_z@description = "Grid Thickness"

    return(del_z)

end




;-------------------------------------------------------------------
; convert from ppm to ug/kg for gas species
; Inputs:
;	species = string of gas species to convert (e.g. "so4")  
;	var_in 	= float/double carrying variable data of any dimensionality

undef("ppm_ugkg_conv")
function ppm_ugkg_conv(species:string,var_in:numeric)
local conv_fac

begin


     ;     should have all gas species details added for completeness...

    
	  if(isStrSubset(species,"_f_c") .or. isStrSubset(species,"_f_o") \
	  		.or. isStrSubset(species,"_b_c") .or. isStrSubset(species,"_b_o")) ; VBS compounds
		   conv_fac = m_dry/m_vbs
	  elseif(isStrSubset(species,"no3"))
		   conv_fac = m_dry/m_no3
	  elseif(isStrSubset(species,"so4"))
		   conv_fac = m_dry/m_so4
	  elseif(isStrSubset(species,"nh4"))
		   conv_fac = m_dry/m_nh4
	  elseif(isStrSubset(species,"cl"))
		   conv_fac = m_dry/m_cl
	  else
	       print("don't recognise species "+species+" for processing in ppm_ugkg_conv")
	       conv_fac = 1.0
	  end if
			   
	  ;; convert from ppm to ppb, then convert to ug/kg		   
	  var_out = (var_in*1e3)/conv_fac
	  return(var_out)
         

end






;
;  function for loading a set of variables from the aerosol phases
;      inputs:
;         var = front string (nh4, so4, na, etc)
;         grp = end string (a01, a02, pm1, pm25, tot, cw) 
;
;  This function assume we're using MOSAIC aerosol, but is suitable for
;  both 4 bin and 8 bin aerosols.
;
undef("load_aero_set_var")
function load_aero_set_var(a:list,var:string,grp:string,time_in:integer)
local temp1, temp2, temp3, temp4, temp5, temp6, temp7, temp8, var_out
begin

	if (isfilevar(a[0],"num_a01"))	; ismosaic

		if (isfilevar(a[0],"num_a08")) ;;; work with 8 bin aerosol
		
			;;; determine which group we're loading
			if(grp .eq. "cw")   ;; loading cloud water variables

				temp1 = (/wrf_user_getvar(a,var+"_cw01",time_in)/)
				temp2 = (/wrf_user_getvar(a,var+"_cw02",time_in)/)
				temp3 = (/wrf_user_getvar(a,var+"_cw03",time_in)/)
				temp4 = (/wrf_user_getvar(a,var+"_cw04",time_in)/)
				temp5 = (/wrf_user_getvar(a,var+"_cw05",time_in)/)
				temp6 = (/wrf_user_getvar(a,var+"_cw06",time_in)/)
				temp7 = (/wrf_user_getvar(a,var+"_cw07",time_in)/)
				temp8 = (/wrf_user_getvar(a,var+"_cw08",time_in)/)
				var_out = temp1 + temp2 + temp3 + temp4 + temp5 + temp6 + temp7 + temp8

			elseif(isStrSubset(grp,"a0")) ;; loading an individual aerosol variable

				var_out = wrf_user_getvar(a,var+grp,time_in)

			elseif( any(grp.eq.(/"pm1","pm25","tot"/)) ) ;;; loading aerosol groups

				temp1 = (/wrf_user_getvar(a,var+"_a01",time_in)/)
				temp2 = (/wrf_user_getvar(a,var+"_a02",time_in)/)
				temp3 = (/wrf_user_getvar(a,var+"_a03",time_in)/)
				temp4 = (/wrf_user_getvar(a,var+"_a04",time_in)/)
				temp5 = (/wrf_user_getvar(a,var+"_a05",time_in)/)

				if(grp .eq. "pm1")  ;; loading pm1 variables  (interstital aerosol only) 
					var_out = temp1 + temp2 + temp3 + temp4 + 0.678*temp5
				elseif(grp .eq. "pm25") ;; loading pm2.5 variables  (interstital aerosol only)
					temp6 = (/wrf_user_getvar(a,var+"_a06",time_in)/)
					var_out = temp1 + temp2 + temp3 + temp4 + temp5 + temp6
				else 
					temp6 = (/wrf_user_getvar(a,var+"_a06",time_in)/)
					temp7 = (/wrf_user_getvar(a,var+"_a07",time_in)/)
					temp8 = (/wrf_user_getvar(a,var+"_a08",time_in)/)
					var_out = temp1 + temp2 + temp3 + temp4 + temp5 + temp6 + temp7 + temp8
				end if	

			else
				print("don't know how to deal with group "+grp+" in load_aero_set_var function")
				exit
			end if

		elseif (isfilevar(a[0],"num_a04")) ;;; work with 4 bin aerosol
		
			;;; determine which group we're loading
			if(grp .eq. "cw")   ;; loading cloud water variables

				temp1 = (/wrf_user_getvar(a,var+"_cw01",time_in)/)
				temp2 = (/wrf_user_getvar(a,var+"_cw02",time_in)/)
				temp3 = (/wrf_user_getvar(a,var+"_cw03",time_in)/)
				temp4 = (/wrf_user_getvar(a,var+"_cw04",time_in)/)
				var_out = temp1 + temp2 + temp3 + temp4

			elseif(isStrSubset(grp,"a0")) ;; loading an individual aerosol variable

				var_out = wrf_user_getvar(a,var+grp,time_in)

			elseif( any(grp.eq.(/"pm1","pm25","tot"/)) ) ;;; loading aerosol groups

				temp1 = (/wrf_user_getvar(a,var+"_a01",time_in)/)
				temp2 = (/wrf_user_getvar(a,var+"_a02",time_in)/)
				temp3 = (/wrf_user_getvar(a,var+"_a03",time_in)/)

				if(grp .eq. "pm1")  ;; loading pm1 variables  (interstital aerosol only) 
					var_out = temp1 + temp2 + (0.678/2.0)*temp3 
				elseif(grp .eq. "pm25") ;; loading pm2.5 variables  (interstital aerosol only)
					var_out = temp1 + temp2 + temp3
				else 
					temp4 = (/wrf_user_getvar(a,var+"_a04",time_in)/)
					var_out = temp1 + temp2 + temp3 + temp4
				end if

			else
				print("don't know how to deal with group "+grp+" in load_aero_set_var function")
				exit
			end if

		
		else
			print("need to use 8 bin or 4 bin MOSAIC aerosol")
			exit
		end if

	else
		print("load_aero_set_var is for MOSAIC aerosol only")
		exit
	end if
	
	
	return(var_out)
	

end  ;;; end load_aero_set_var





;
;   Function for calculating the aerosol pH
;
;
undef("get_aerosol_ph")
function get_aerosol_ph(a:list,grp:string,time_in:integer)
local so4_mol, no3_mol, nh4_mol, na_mol, cl_mol, hplu_mol, aero_mass, h_molarity, ph, \
	so4_mass, no3_mass, nh4_mass, na_mass, cl_mass, h2o_mass, \
	so4_vol, no3_vol, nh4_vol, na_vol, cl_vol, h2o_vol
begin

	;;; load chemical species & convert to moles / kg_dry_air (with conversion from ug to g for the mass)
	so4_mass = load_aero_set_var(a,"so4",grp,time_in) *1e-6
	so4_mol  = so4_mass  / m_so4
	no3_mass = load_aero_set_var(a,"no3",grp,time_in) *1e-6
	no3_mol  = no3_mass  / m_no3 
	nh4_mass = load_aero_set_var(a,"nh4",grp,time_in) *1e-6
	nh4_mol  = nh4_mass  / m_nh4 
	na_mass  = load_aero_set_var(a,"na",grp,time_in) *1e-6
	na_mol   = na_mass  / m_na 
	cl_mass  = load_aero_set_var(a,"cl",grp,time_in) *1e-6
	cl_mol   = cl_mass  / m_cl

	h2o_mass = load_aero_set_var(a,"water",grp,time_in) *1e-6

	;;; calculate the number of moles of H+ 
	hplu_mol = 2.0*so4_mol + no3_mol + cl_mol - nh4_mol - na_mol	
	;hplu_mol = where(hplu_mol .gt. 0.0, hplu_mol, 0.0)

	;;; calculate the volumes of each aerosol component (in cm3)
	so4_vol = so4_mass / dens_so4_aer
	no3_vol = no3_mass / dens_no3_aer
	nh4_vol = nh4_mass / dens_nh4_aer
	na_vol  = na_mass  / dens_na_aer
	cl_vol  = cl_mass  / dens_cl_aer
	h2o_vol = h2o_mass / dens_water_aer
	
	;;; determine total volume & convert to litres (assume H+ contribution to mass is insignificant)
	total_volume = (so4_vol+no3_vol+nh4_vol+na_vol+cl_vol+h2o_vol) *1e-3

	
	;;; calculate H+ molarity
	h_molarity = ( hplu_mol + 1e-7*(h2o_vol*1e-3)) / total_volume
	
	hplu_stats = calc_standard_stats(h_molarity)
	
	print(hplu_stats)

	
	;;; if the h_molarity = 0 then set it to -1, so we have NaN values for pH, not infinity
	h_molarity = where(h_molarity .gt. 0.0, h_molarity, 0.0)
	
	;;; calculate pH
	ph = new(dimsizes(h_molarity),float,-9999)
	ph = where(h_molarity .gt. 0.0 , -log10(h_molarity), -9999)
	
	;;; transfer file atts
	if(grp.eq."pm1")
		ph@description = "PM1 pH"
	elseif(grp.eq."pm25")
		ph@description = "PM2.5 pH"
	elseif(grp.eq."tot")
		ph@description = "PM10 pH"
	else
		ph@description = "aerosol "+grp+" pH"
	end if
	ph@units = ""


	return(ph)

end



;
;   Function for calculating cloud water pH
;
;   This calculates the number of H+ moles (using ion balance for cw variables)
;     and then uses the QCLOUD variable to get the liquid volume for calculating
;     H+ molarity  (pH = -log10(H+ molarity) )
;
undef("get_cloud_ph")
function get_cloud_ph(a:list,time_in:integer)
local so4_mol, no3_mol, nh4_mol, na_mol, cl_mol, hplu_mol, qcloud, h_molarity, ph
begin

	;;; load chemical species & convert to moles / kg_dry_air (with conversion from ug to g)
	so4_mol = load_aero_set_var(a,"so4","cw",time_in) *1e-6 / m_so4
	no3_mol = load_aero_set_var(a,"no3","cw",time_in) *1e-6 / m_no3 
	nh4_mol = load_aero_set_var(a,"nh4","cw",time_in) *1e-6 / m_nh4 
	na_mol = load_aero_set_var(a,"na","cw",time_in) *1e-6 / m_na 
	cl_mol = load_aero_set_var(a,"cl","cw",time_in) *1e-6 / m_cl

	;;; calculate the number of moles of H+ 
	hplu_mol = 2.0*so4_mol + no3_mol + cl_mol - nh4_mol - na_mol
	;hplu_mol = where(hplu_mol .gt. 0.0, hplu_mol, 0.0)

	;;; load the cloud water loading (and set zero to -1, to avoid division errors)
	qcloud = (/wrf_user_getvar(a,"QCLOUD",time_in)/)
	qcloud = where(qcloud .gt. 0.0, qcloud, -1.0)

	;;; calculate H+ molarity
	h_molarity = new(dimsizes(qcloud),float,-9999)
	h_molarity = where(qcloud .gt. 0.0, ( hplu_mol + 1e-7*(qcloud)) / qcloud , -9999.)

	hplu_stats = calc_standard_stats(h_molarity)
	
	print(hplu_stats)

	
	;;; calculate pH
	ph = new(dimsizes(h_molarity),float,-9999)
	ph = where(h_molarity .gt. 0.0 , -log10(h_molarity), -9999)
	
	;;; transfer file atts
	ph@description = "Cloud Water pH"
	ph@units = ""


	return(ph)

end






;
;   Function for calculating approximate rain water pH
;
;   This calculates the number of H+ moles (using ion balance for cw variables)
;     and then uses the QCLOUD variable to get the liquid volume for calculating
;     H+ molarity  (pH = -log10(H+ molarity) )
;
;   Then we obtain a weighted average of the Molarity above anywhere which has had
;     rain in the last hour (this is an instantaneous average, so only roughly 
;     applicable to the rain itself).
;
;
undef("get_rain_approx_ph")
function get_rain_approx_ph(a:list,time_in:integer)
local time, so4_mol, no3_mol, nh4_mol, na_mol, cl_mol, 
		hplu_mol, qcloud, data_dim, \
		h_molarity_cloud, PH, PHB, height, zdiff, kk, \
		h_mol_weight_av, weight_factor, cloud_weight, tk, ik, jk, \
		tot_w_fac, var1, var2, var3, rain_accum, rain_hour, \
		ph 
begin

	time = time_in(0)  ; get some time information

	;;; load chemical species & convert to moles / kg_dry_air (with conversion from ug to g)
	so4_mol = load_aero_set_var(a,"so4","cw",time_in) *1e-6 / m_so4
	no3_mol = load_aero_set_var(a,"no3","cw",time_in) *1e-6 / m_no3 
	nh4_mol = load_aero_set_var(a,"nh4","cw",time_in) *1e-6 / m_nh4 
	na_mol = load_aero_set_var(a,"na","cw",time_in) *1e-6 / m_na 
	cl_mol = load_aero_set_var(a,"cl","cw",time_in) *1e-6 / m_cl

	;;; calculate the number of moles of H+ 
	hplu_mol = 2.0*so4_mol + no3_mol + cl_mol - nh4_mol - na_mol
	;hplu_mol = where(hplu_mol .gt. 0.0, hplu_mol, 0.0)

	;;; load the cloud water loading (and set zero to -1, to avoid division errors)
	qcloud = (/wrf_user_getvar(a,"QCLOUD",time_in)/)
	qcloud = where(qcloud .gt. 0.0, qcloud, -1.0)

	;;; get the dimensions of qcloud
	data_dim = dimsizes(qcloud)

	;;; calculate H+ molarity (in whole volume of atmosphere)
	h_molarity_cloud = new(data_dim,float,-9999)
	h_molarity_cloud = where(qcloud .gt. 0.0, ( hplu_mol + 1e-7*(qcloud)) / qcloud , -9999.)


	;;; obtain the layer thicknesses - calculating the heights at layer boundaries
	if(time.eq.-1)
		PH  = a[:]->PH     ; staggered on Z 
		PHB = a[:]->PHB    ; staggered on Z
	else
		PH  = a[:]->PH(time_in,:,:,:)     ; staggered on Z 
		PHB = a[:]->PHB(time_in,:,:,:)    ; staggered on Z	
	end if
	height = (PH + PHB)/9.8
	zdiff  = new(data_dim,float,-9999) ; unstaggered dataset
	do kk = 0,data_dim(1)-1
		zdiff(:,kk,:,:) = (height(:,kk+1,:,:) - height(:,kk,:,:))
	end do
	delete([/PH,PHB,height/])

	
	;;; calculate the weight for each cloudy layer in a column
	weight_factor = where(qcloud .gt. 0.0, qcloud*zdiff, 0.0)  ; this is dimensionless
	cloud_weight = new(data_dim,float,-9999)
	do tk = 0,data_dim(0)-1
	do ik = 0,data_dim(2)-1
	do jk = 0,data_dim(3)-1
		tot_w_fac = sum(weight_factor(tk,:,ik,jk))
		if(tot_w_fac.gt.0.0)
			cloud_weight(tk,:,ik,jk) = weight_factor(tk,:,ik,jk)/tot_w_fac
		else
			cloud_weight(tk,:,ik,jk) = 0.0
		end if
	end do
	end do
	end do
	
	;;; calculate a weighted H+ molarity through the model column
	h_mol_weight_av = new((/data_dim(0),data_dim(2),data_dim(3)/),float,-9999)	
	do tk = 0,data_dim(0)-1
	do ik = 0,data_dim(2)-1
	do jk = 0,data_dim(3)-1
		h_mol_weight_av(tk,ik,jk) = sum(h_molarity_cloud(tk,:,ik,jk)*cloud_weight(tk,:,ik,jk))
	end do
	end do
	end do

	
	;;; get the accumulated rain fall
	var1 = a[:]->RAINC	
	var2 = a[:]->RAINSH	
	var3 = a[:]->RAINNC	
	rain_accum = var1+var2+var3
	
	;;; determine the hourly rainfall
	rain_hour = rain_accum * 0.0
	do tk= 1,data_dim(0)-1
			rain_hour(tk,:,:) = rain_accum(tk,:,:) - rain_accum(tk-1,:,:)
	end do
	
	;;; calculate pH
	ph = new(dimsizes(h_mol_weight_av),float,-9999)
	ph = where(h_mol_weight_av .gt. 0.0 , -log10(h_mol_weight_av), -9999)
	
	;; filter pH data, to match rainfall locations
	ph = where(rain_hour.gt.0.0,ph,-9999)

	;;; transfer file atts
	ph@description = "Estimated rainfall pH"
	ph@units = ""

	return(ph)

end












;-------------------------------------------------------------------;
; AEROSOL SUMMING FUNCTIONS
; SAN - 17/06/14: scripts modified to not process meta data, as this improves efficiency.
;                 syntax: x = (/ y /) -> x will not have the meta_data of y.

;-------------------------------------------------------------------;
;     get_aerosol     -     extracts and sums up dry aerosol species into one variable

; Using following definitions for aerosol sizes (in 8-bin mode):
; fine 	= smallest 4 bins ( Dp < 0.675um)
; coarse= largest 2 bins (Dp > 2.5 um)
; tot 	= All bins
; pm25	= Smallest 6 bins (Dp < 2.5um)
; pm1	= smallest 4 and a bit bins (dp < 1um)

; SAN - 23/01/14: changed definition of coarse mode to >2.5um.

; SAN - 20/04/15: adding options for use with GOCART aerosol mechanism
;  Can only be used to find the total aerosol and sulfate in GOCART mode.

; SAN - 21/04/15: redesigned function to make neater and expand functionality:
;  * Merged with get_aerosol_full, set inc_cw = True to use in same way as get_aerosol_full
;  * New argument inc_cw: set to True to include CW variables, otherwise set to False.
;						 Does nothing if running with GOCART aerosol.
;  * Makes arrays containing which bins to sum, then loops over all bins. 
;  * Can work with diagnostic MOSAIC VBS species.

undef("get_aerosol")
function get_aerosol(a:list,species:string,inc_cw:logical,time_in:integer)
local a1, sh_spec, spec_list, nspec, dummy, temp, i, desc, aero_mech, nbins, nphase, \
	fine_array, coarse_array, tot_array, pm25_array, pm1_array, array_map, mode_name, \
	bin_names, aq_string, conv_fac
begin

  ; Do quick and dirty check to see if input species name is an aerosol species.
  if(isStrSubset(species,"_") )   

	; check if the mechanism is GOCART or MOSAIC
	if (isfilevar(a[0],"SEAS_1"))	; isgocart
	
;**********************     GOCART aerosol     *************************;

		aero_mech = "GOCART" 
		print("get_aerosol: summing GOCART species")
		; Declare conv_fac, this is a dummy variable except for with sulf.
		conv_fac = 1.
		
		; Create a list of species to sum
			
		; OC and BC - sum hydrophobic and hydrophilic components
		if ( any(species .eq. (/"oc_tot", "bc_tot" /)) )
		
        	sh_spec = str_sub_str(species, "c_tot","c" )
        	sh_spec = str_upper(sh_spec)   	
        	spec_list = (/sh_spec+"1",sh_spec+"2"/)
        	desc = "Total "+sh_spec
        
		elseif ( any(species .eq. (/ "seas_tot", "SEAS_tot"/)) )
			; 4 sea-salt bins
        	spec_list = (/"SEAS_1","SEAS_2","SEAS_3","SEAS_4"/)
			desc      = "Total sea-salt"

		elseif ( any(species .eq. (/ "dust_tot", "DUST_tot"/)) )
			; 5 dust bins
        	spec_list = (/"DUST_1","DUST_2","DUST_3","DUST_4", "DUST_5"/)
			desc      = "Total dust"
		

		elseif(species .eq. "sulf_tot")
			; estimate the total sulfate mass, including ammonium componant
			; See sum_pm_gocart in module_gocart_aerosols.F for equation		
			spec_list = (/ "sulf" /)
			desc = "Total Ammonium sulfate"
			
			; Define conversion factor for ppm->ug/kg, include estimate for NH4 componant  
			conv_fac = 1e3 * (m_so4/m_dry) *(2*m_nh4 + m_s04)/m_so4		
			
		else 
		
			print("get_aerosol: can only sum OC_tot, BC_tot, dust_tot, seas_tot and sulf_tot GOCART aerosol")
			break
			
		end if
		
		nspec = dimsizes(spec_list)
		var_a = wrf_user_getvar(a,spec_list(0),time_in) * conv_fac
		
		if (nspec .gt. 0)
			do i = 1, nspec-1
				var_a = var_a + (/ wrf_user_getvar(a, spec_list(i), time_in) /) * conv_fac
			end do
		end if
		
		var_a@units       = "ug/kg-dryair"
		var_a@description = desc
		
;**********************     MOSAIC aerosol     *************************;
	
	elseif (isfilevar(a[0],"num_a01"))	; ismosaic

		; Declare mapping arrays (in 8-bin mode). 
		fine_array   = (/ 1., 1., 1., 1., 0., 0., 0., 0. /)	; fine = smallest 4 bins
		coarse_array = (/ 0., 0., 0., 0., 0., 0., 1., 1. /) ; coarse = largest 2 bins
		tot_array    = (/ 1., 1., 1., 1., 1., 1., 1., 1. /) ; tot  = all bins
		pm25_array   = (/ 1., 1., 1., 1., 1., 1., 0., 0. /) ; pm25 = smallest 6 bins
		pm1_array    = (/ 1., 1., 1., 1., 0.678, 0., 0., 0. /) ; pm1 = smallest 4 bins + 68% of bin 5
		a01_array    = (/ 1., 0., 0., 0., 0., 0., 0., 0. /)
		a02_array    = (/ 0., 1., 0., 0., 0., 0., 0., 0. /)
		a03_array    = (/ 0., 0., 1., 0., 0., 0., 0., 0. /)
		a04_array    = (/ 0., 0., 0., 1., 0., 0., 0., 0. /)
		a05_array    = (/ 0., 0., 0., 0., 1., 0., 0., 0. /)
		a06_array    = (/ 0., 0., 0., 0., 0., 1., 0., 0. /)
		a07_array    = (/ 0., 0., 0., 0., 0., 0., 1., 0. /)
		a08_array    = (/ 0., 0., 0., 0., 0., 0., 0., 1. /)
		cw01_array    = (/ 1., 0., 0., 0., 0., 0., 0., 0. /)
		cw02_array    = (/ 0., 1., 0., 0., 0., 0., 0., 0. /)
		cw03_array    = (/ 0., 0., 1., 0., 0., 0., 0., 0. /)
		cw04_array    = (/ 0., 0., 0., 1., 0., 0., 0., 0. /)
		cw05_array    = (/ 0., 0., 0., 0., 1., 0., 0., 0. /)
		cw06_array    = (/ 0., 0., 0., 0., 0., 1., 0., 0. /)
		cw07_array    = (/ 0., 0., 0., 0., 0., 0., 1., 0. /)
		cw08_array    = (/ 0., 0., 0., 0., 0., 0., 0., 1. /)
		
		if(isStrSubset(species,"fine") ) ; returns true if species contais "fine"		 
			; refine species name
			sh_spec   = str_sub_str(species, "_fine","" )
			array_map = fine_array
			mode_name = "fine mode"
		elseif(isStrSubset(species,"coarse") ) ; returns true if species contais "coarse"
			sh_spec   = str_sub_str(species, "_coarse","" )
			array_map = coarse_array
			mode_name = "coarse mode"
		elseif(isStrSubset(species,"tot") ) ; returns true if species contais "tot"		 
			sh_spec = str_sub_str(species, "_tot","" )
			array_map = tot_array
			mode_name = "PM10"
		elseif(isStrSubset(species,"pm25") ) ; returns true if species contais "pm25"		 
			sh_spec = str_sub_str(species, "_pm25","" )
			array_map = pm25_array
			mode_name = "PM2.5"
		elseif(isStrSubset(species,"pm1") ) ; returns true if species contais "pm1"
			sh_spec   = str_sub_str(species, "_pm1","" )
			array_map = pm1_array
			mode_name = "PM1"
		elseif(isStrSubset(species,"a01") )
			sh_spec   = str_sub_str(species, "_a01","" )
			array_map = a01_array
			mode_name = "a01"
		elseif(isStrSubset(species,"a02") )
			sh_spec   = str_sub_str(species, "_a02","" )
			array_map = a02_array
			mode_name = "a02"
		elseif(isStrSubset(species,"a01") )
			sh_spec   = str_sub_str(species, "_a01","" )
			array_map = a01_array
			mode_name = "a01"
		elseif(isStrSubset(species,"a03") )
			sh_spec   = str_sub_str(species, "_a03","" )
			array_map = a03_array
			mode_name = "a03"
		elseif(isStrSubset(species,"a04") )
			sh_spec   = str_sub_str(species, "_a04","" )
			array_map = a04_array
			mode_name = "a04"
		elseif(isStrSubset(species,"a05") )
			sh_spec   = str_sub_str(species, "_a05","" )
			array_map = a05_array
			mode_name = "a05"
		elseif(isStrSubset(species,"a06") )
			sh_spec   = str_sub_str(species, "_a06","" )
			array_map = a06_array
			mode_name = "a06"
		elseif(isStrSubset(species,"a07") )
			sh_spec   = str_sub_str(species, "_a07","" )
			array_map = a07_array
			mode_name = "a07"
		elseif(isStrSubset(species,"a08") )
			sh_spec   = str_sub_str(species, "_a08","" )
			array_map = a08_array
			mode_name = "a08"
		elseif(isStrSubset(species,"cw01") )
			sh_spec   = str_sub_str(species, "_cw01","" )
			array_map = cw01_array
			mode_name = "cw01"
		elseif(isStrSubset(species,"cw02") )
			sh_spec   = str_sub_str(species, "_cw02","" )
			array_map = cw02_array
			mode_name = "cw02"
		elseif(isStrSubset(species,"cw01") )
			sh_spec   = str_sub_str(species, "_cw01","" )
			array_map = cw01_array
			mode_name = "cw01"
		elseif(isStrSubset(species,"cw03") )
			sh_spec   = str_sub_str(species, "_cw03","" )
			array_map = cw03_array
			mode_name = "cw03"
		elseif(isStrSubset(species,"cw04") )
			sh_spec   = str_sub_str(species, "_cw04","" )
			array_map = cw04_array
			mode_name = "cw04"
		elseif(isStrSubset(species,"cw05") )
			sh_spec   = str_sub_str(species, "_cw05","" )
			array_map = cw05_array
			mode_name = "cw05"
		elseif(isStrSubset(species,"cw06") )
			sh_spec   = str_sub_str(species, "_cw06","" )
			array_map = cw06_array
			mode_name = "cw06"
		elseif(isStrSubset(species,"cw07") )
			sh_spec   = str_sub_str(species, "_cw07","" )
			array_map = cw07_array
			mode_name = "cw07"
		elseif(isStrSubset(species,"cw08") )
			sh_spec   = str_sub_str(species, "_cw08","" )
			array_map = cw08_array
			mode_name = "cw08"
		else
			print("get_aerosol: ERROR - input name must contain 'fine', 'coarse', 'tot', 'pm25' or 'pm1'.")
			exit			
		end if
			
		; MOSAIC VBS diagnostic species use UC "A" instead of lc "a".
		if (any(sh_spec .eq. (/"HOA", "SOA", "BBOA", "BBSOA", "HSOA", "BIOG", "TOTOA"/)) )
			; note no AQ variables for VBS diagnostic species
			bin_names = (/ (/"A01", "A02", "A03", "A04", "A05", "A06", "A07", "A08" /), \
			               (/"N/A", "N/A", "N/A", "N/A", "N/A", "N/A", "N/A", "N/A" /)  /)
		; Special case: summed GLYSOA species have              
		elseif(sh_spec .eq. "GLYSOA")
			bin_names = (/ (/"A01", "A02", "A03", "A04", "A05", "A06", "A07", "A08" /), \
			               (/"CW01","CW02","CW03","CW04","CW05","CW06","CW07","CW08"/)  /)
		else
			bin_names = (/ (/"a01", "a02", "a03", "a04", "a05", "a06", "a07", "a08" /), \
			               (/"cw01","cw02","cw03","cw04","cw05","cw06","cw07","cw08"/)  /)
		end if
		
		; check if there is Aqueous chemistry or not
		if (inc_cw)
			nphase = 2
			aq_string = " including aqueous phase."
			
			; check if there are cloud water variables
			if (.not. isfilevar(a[0],"num_cw01"))
				print("get_aerosol: ERROR - Cannot sum CW variables")
				print("Either run with aqueous chemistry or set inc_cw = False")
				exit			
			end if
			
			if (bin_names(1,0) .eq. "N/A") 
				print("get_aerosol: no CW variables for "+sh_spec)
				print("Cloud water variables already included in diagnostic species.")
				nphase = 1
			end if
		else
			nphase = 1
			aq_string = "."
		end if
		
		; Load a dummy for copying attributes
		dummy = wrf_user_getvar(a,sh_spec+"_"+bin_names(0,0),time_in)
				
		; Check if there is an 8th bin, if so add in 8 bin mode, else 4-bin mode
		if (isfilevar(a[0],"num_a08")) ; is8bin
			aero_mech = "MOSAIC_8bin"
			nbin = 8		
			
		; Otherwise sum up in 4-bin mode
		elseif (isfilevar(a[0],"num_a04")) ; is4bin
			aero_mech = "MOSAIC_4bin"
			nbin = 4

			; resize the array map
			temp = new(4, float)
			do i = 0, 3
				temp(i) = avg(array_map(2*i:2*i+1))
			end do
			delete(array_map)
			array_map = temp
			delete(temp)
			
		else
			print("get_aerosol: ERROR - Cannot find bin 4, currently written for MOSAIC 4 or 8bin")
			exit
		end if
		
		print("get_aerosol: summing "+species+" in "+aero_mech+" mode") 
		; print(array_map)
		
		var_a = (/ wrf_user_getvar(a,sh_spec+"_"+bin_names(0,0),time_in) /) * 0.0
		
		; Loop across all size bins, summing the aerosol
		do iphase = 0, nphase-1
			do ibin = 0, nbin-1
				; Only load data from file if nonzero
				if (array_map(ibin) .gt. 0.0)
					if( isfilevar(a[0],sh_spec+"_"+bin_names(iphase,ibin)))
						a1 = (/ wrf_user_getvar(a,sh_spec+"_"+bin_names(iphase,ibin),time_in) /)
						var_a = var_a + a1*array_map(ibin)
						delete(a1)
					else
						print(sh_spec+"_"+bin_names(iphase,ibin)+" does not exist, skipping...")
					end if
				else
					continue
				end if	
			end do
		end do
				
		; Copy meta data from dummy	
		copy_VarMeta (dummy, var_a) 

		; Modify description:
		desc = dummy@description
		if (isStrSubset(desc, "bin 01"))
			desc = str_sub_str(desc, "bin 01", mode_name)
		else
			desc = sh_spec +" "+mode_name
		end if
		
		var_a@description = desc + aq_string   
		         
	else
		print("get_aerosol: ERROR - Currently only written for GOCART and MOSAIC 4/8bin aerosol mechanisms")
		exit
	end if

  else
	  print("get_aerosol: "+species+" is not an aerosol species, needs '_' in name")
	  break
  end if

  return(var_a)

end


;----------------------------------------------------------------------;
;
;  get_vbs_aerosol - extracts and sums up the VBS species into one variable
;
;      This uses the get_aerosol function to extract data
; DL 8/5/2015 - using inc_cw to control inclusion of cloud variables. Also added more options
;               for VBS schemes to extract (and a check for what variables are in the file).

undef("vbs_species_extract_loop")
function vbs_species_extract_loop(a:list,species:string,vbs_string:string,vbs_species:string,\
								inc_cw:logical,time_in:integer)
local var_a, num_vbs, ispec, spc_extract, dummy
begin

	; record number of vbs species to be extracted
	num_vbs = dimsizes(vbs_species)

	; loop through the species listed for extraction
	spc_extract = str_sub_str(species, vbs_string, vbs_species(0))
	var_a = (/ get_aerosol(a,spc_extract,inc_cw,time_in) /)
	do ispec = 1, num_vbs-1
		spc_extract = str_sub_str(species, vbs_string, vbs_species(ispec))
		dummy = (/ get_aerosol(a,spc_extract,inc_cw,time_in) /)
		var_a = var_a + dummy
	end do


	return(var_a)

end



undef("get_vbs_aerosol")
function get_vbs_aerosol(a:list,species:string,inc_cw:logical,time_in:integer)
local vbs_species, vbs_species_exist, vbs_string, num_vbs, ispec, dummy, units, var_b, desc, var_a, var_c
begin


    if(isStrSubset(species,"_") )     ; only have "_" within names of aerosol species
		
		units = "ug/kg-dryair"
		
		if(isStrSubset(species,"amsratio")) ; solve inorganic:organic ratio

			vbs_string = "amsratio"
			
			;;;;; first extract the organic data
			vbs_species = (/ "pcg1_f_c", "pcg2_f_c", "pcg3_f_c", "pcg4_f_c", \
							 "pcg5_f_c", "pcg6_f_c", "pcg7_f_c", "pcg8_f_c", "pcg9_f_c", \
							 "pcg1_b_c", "pcg2_b_c", "pcg3_b_c", "pcg4_b_c", \
							 "pcg5_b_c", "pcg6_b_c", "pcg7_b_c", "pcg8_b_c", "pcg9_b_c", \
							 "opcg1_f_c", "opcg2_f_c", "opcg3_f_c", "opcg4_f_c", \
							 "opcg5_f_c", "opcg6_f_c", "opcg7_f_c", "opcg8_f_c", \
							 "opcg1_b_c", "opcg2_b_c", "opcg3_b_c", "opcg4_b_c", \
							 "opcg5_b_c", "opcg6_b_c", "opcg7_b_c", "opcg8_b_c", \
							 "pcg1_f_o", "pcg2_f_o", "pcg3_f_o", "pcg4_f_o", \
							 "pcg5_f_o", "pcg6_f_o", "pcg7_f_o", "pcg8_f_o", "pcg9_f_o", \
							 "pcg1_b_o", "pcg2_b_o", "pcg3_b_o", "pcg4_b_o", \
							 "pcg5_b_o", "pcg6_b_o", "pcg7_b_o", "pcg8_b_o", "pcg9_b_o", \
							 "opcg1_f_o", "opcg2_f_o", "opcg3_f_o", "opcg4_f_o", \
							 "opcg5_f_o", "opcg6_f_o", "opcg7_f_o", "opcg8_f_o", \
							 "opcg1_b_o", "opcg2_b_o", "opcg3_b_o", "opcg4_b_o", \
							 "opcg5_b_o", "opcg6_b_o", "opcg7_b_o", "opcg8_b_o", \
							 "oc" /)
			; obtain the organic mass
			var_a = vbs_species_extract_loop(a,species,vbs_string,vbs_species,inc_cw,time_in)


			;;;;; then extract the inorganic data
			vbs_species := (/ "so4", "no3", "nh4", "na", "cl" /)
			; obtain the inorganic mass
			var_b = vbs_species_extract_loop(a,species,vbs_string,vbs_species,inc_cw,time_in)


			;;;;; now calculate the inorganic:organic ratio
			
			var_a = var_b/where(var_a .gt. 1e-12, var_a, 1e-12)


			units = ""

		elseif(isStrSubset(species,"ioratio")) ; solve inorganic:organic ratio

			vbs_string = "ioratio"
			
			;;;;; first extract the organic data
			vbs_species = (/ "pcg1_f_c", "pcg2_f_c", "pcg3_f_c", "pcg4_f_c", \
							 "pcg5_f_c", "pcg6_f_c", "pcg7_f_c", "pcg8_f_c", "pcg9_f_c", \
							 "pcg1_b_c", "pcg2_b_c", "pcg3_b_c", "pcg4_b_c", \
							 "pcg5_b_c", "pcg6_b_c", "pcg7_b_c", "pcg8_b_c", "pcg9_b_c", \
							 "opcg1_f_c", "opcg2_f_c", "opcg3_f_c", "opcg4_f_c", \
							 "opcg5_f_c", "opcg6_f_c", "opcg7_f_c", "opcg8_f_c", \
							 "opcg1_b_c", "opcg2_b_c", "opcg3_b_c", "opcg4_b_c", \
							 "opcg5_b_c", "opcg6_b_c", "opcg7_b_c", "opcg8_b_c", \
							 "pcg1_f_o", "pcg2_f_o", "pcg3_f_o", "pcg4_f_o", \
							 "pcg5_f_o", "pcg6_f_o", "pcg7_f_o", "pcg8_f_o", "pcg9_f_o", \
							 "pcg1_b_o", "pcg2_b_o", "pcg3_b_o", "pcg4_b_o", \
							 "pcg5_b_o", "pcg6_b_o", "pcg7_b_o", "pcg8_b_o", "pcg9_b_o", \
							 "opcg1_f_o", "opcg2_f_o", "opcg3_f_o", "opcg4_f_o", \
							 "opcg5_f_o", "opcg6_f_o", "opcg7_f_o", "opcg8_f_o", \
							 "opcg1_b_o", "opcg2_b_o", "opcg3_b_o", "opcg4_b_o", \
							 "opcg5_b_o", "opcg6_b_o", "opcg7_b_o", "opcg8_b_o", \
							 "oc" /)
			; obtain the organic mass
			var_a = vbs_species_extract_loop(a,species,vbs_string,vbs_species,inc_cw,time_in)


			;;;;; then extract the inorganic data
			vbs_species := (/ "bc", "so4", "no3", "nh4", "na", "cl", "oin" /)
			; obtain the inorganic mass
			var_b = vbs_species_extract_loop(a,species,vbs_string,vbs_species,inc_cw,time_in)

			;;;;; now calculate the inorganic:organic ratio
			
			var_a = var_b/where(var_a .gt. 1e-12, var_a, 1e-12)


			units = ""

		elseif(isStrSubset(species,"vbs_soa_ocratio")) ; solve O:C ratio

			vbs_string = "vbs_soa_ocratio"
			
			;;;;; first extract the Carbon data
			vbs_species = (/ "pcg1_f_c", "pcg2_f_c", "pcg3_f_c", "pcg4_f_c", \
							 "pcg5_f_c", "pcg6_f_c", "pcg7_f_c", "pcg8_f_c", "pcg9_f_c", \
							 "pcg1_b_c", "pcg2_b_c", "pcg3_b_c", "pcg4_b_c", \
							 "pcg5_b_c", "pcg6_b_c", "pcg7_b_c", "pcg8_b_c", "pcg9_b_c", \
							 "opcg1_f_c", "opcg2_f_c", "opcg3_f_c", "opcg4_f_c", \
							 "opcg5_f_c", "opcg6_f_c", "opcg7_f_c", "opcg8_f_c", \
							 "opcg1_b_c", "opcg2_b_c", "opcg3_b_c", "opcg4_b_c", \
							 "opcg5_b_c", "opcg6_b_c", "opcg7_b_c", "opcg8_b_c" \
							/)
			; obtain the carbon mass
			var_a = vbs_species_extract_loop(a,species,vbs_string,vbs_species,inc_cw,time_in)

			; convert from "carbon" to proper carbon (assuming a fixed carbon:other ratio in model)
			;        this is taken directly from module_mosaic_driver.F
			var_a = var_a * (180.0/211.0)


			;;;;; then extract the Oxygen data
			vbs_species = (/ "pcg1_f_o", "pcg2_f_o", "pcg3_f_o", "pcg4_f_o", \
							 "pcg5_f_o", "pcg6_f_o", "pcg7_f_o", "pcg8_f_o", "pcg9_f_o", \
							 "pcg1_b_o", "pcg2_b_o", "pcg3_b_o", "pcg4_b_o", \
							 "pcg5_b_o", "pcg6_b_o", "pcg7_b_o", "pcg8_b_o", "pcg9_b_o", \
							 "opcg1_f_o", "opcg2_f_o", "opcg3_f_o", "opcg4_f_o", \
							 "opcg5_f_o", "opcg6_f_o", "opcg7_f_o", "opcg8_f_o", \
							 "opcg1_b_o", "opcg2_b_o", "opcg3_b_o", "opcg4_b_o", \
							 "opcg5_b_o", "opcg6_b_o", "opcg7_b_o", "opcg8_b_o" \
							/)
			; obtain the oxygen mass
			var_b = vbs_species_extract_loop(a,species,vbs_string,vbs_species,inc_cw,time_in)

			;;;;; then extract the SOA data
			vbs_species := (/ "ant1_c", "ant2_c", "ant3_c", "ant4_c", \
							 "biog1_c", "biog2_c", "biog3_c", "biog4_c" /)
			; obtain the SOA mass
			var_c = vbs_species_extract_loop(a,species,vbs_string,vbs_species,inc_cw,time_in)

			; add carbon from SOA to carbon data (converting from SOA mass to carbon mass, following Shrivastava et al, 2011)
			var_a = var_a + (/ var_c / 1.9 /)  ; 1.9 = ratio of OM to OC
			; add oxygen from SOA to oxygen data (converting from SOA mass to oxygen mass, following Shrivastava et al, 2011)
			var_b = var_b + (/ var_c / (1.9*1.875) /)  ; 1.9 = ratio of OM to OC, and 1.875 = C:O mass ratio


			;;;;; now calculate the O:C ratio
			; convert from ug/kg_air to ~mol/kg_air
			var_a = var_a / m_carbon
			var_b = var_b / m_oxygen
			
			var_a = var_b/where(var_a .gt. 1e-12, var_a, 1e-12)


			units = ""

		elseif(isStrSubset(species,"anthro_ocratio")) ; solve O:C ratio

			vbs_string = "anthro_ocratio"
			
			;;;;; first extract the Carbon data
			vbs_species = (/ "pcg1_f_c", "pcg2_f_c", "pcg3_f_c", "pcg4_f_c", \
							 "pcg5_f_c", "pcg6_f_c", "pcg7_f_c", "pcg8_f_c", "pcg9_f_c", \
							 "opcg1_f_c", "opcg2_f_c", "opcg3_f_c", "opcg4_f_c", \
							 "opcg5_f_c", "opcg6_f_c", "opcg7_f_c", "opcg8_f_c" \
							/)
			; obtain the carbon mass
			var_a = vbs_species_extract_loop(a,species,vbs_string,vbs_species,inc_cw,time_in)

			; convert from "carbon" to proper carbon (assuming a fixed carbon:other ratio in model)
			;        this is taken directly from module_mosaic_driver.F
			var_a = var_a * (180.0/211.0)


			;;;;; then extract the Oxygen data
			vbs_species = (/ "pcg1_f_o", "pcg2_f_o", "pcg3_f_o", "pcg4_f_o", \
							 "pcg5_f_o", "pcg6_f_o", "pcg7_f_o", "pcg8_f_o", "pcg9_f_o", \
							 "opcg1_f_o", "opcg2_f_o", "opcg3_f_o", "opcg4_f_o", \
							 "opcg5_f_o", "opcg6_f_o", "opcg7_f_o", "opcg8_f_o" \
							/)
			; obtain the oxygen mass
			var_b = vbs_species_extract_loop(a,species,vbs_string,vbs_species,inc_cw,time_in)

			;;;;; now calculate the O:C ratio
			; convert from ug/kg_air to ~mol/kg_air
			var_a = var_a / m_carbon
			var_b = var_b / m_oxygen
			
			var_a = var_b/where(var_a .gt. 1e-12, var_a, 1e-12)


			units = ""

		elseif(isStrSubset(species,"ocratio")) ; solve O:C ratio

			vbs_string = "ocratio"
			
			;;;;; first extract the Carbon data
			vbs_species = (/ "pcg1_f_c", "pcg2_f_c", "pcg3_f_c", "pcg4_f_c", \
							 "pcg5_f_c", "pcg6_f_c", "pcg7_f_c", "pcg8_f_c", "pcg9_f_c", \
							 "pcg1_b_c", "pcg2_b_c", "pcg3_b_c", "pcg4_b_c", \
							 "pcg5_b_c", "pcg6_b_c", "pcg7_b_c", "pcg8_b_c", "pcg9_b_c", \
							 "opcg1_f_c", "opcg2_f_c", "opcg3_f_c", "opcg4_f_c", \
							 "opcg5_f_c", "opcg6_f_c", "opcg7_f_c", "opcg8_f_c", \
							 "opcg1_b_c", "opcg2_b_c", "opcg3_b_c", "opcg4_b_c", \
							 "opcg5_b_c", "opcg6_b_c", "opcg7_b_c", "opcg8_b_c" \
							/)
			; obtain the carbon mass
			var_a = vbs_species_extract_loop(a,species,vbs_string,vbs_species,inc_cw,time_in)

			; convert from "carbon" to proper carbon (assuming a fixed carbon:other ratio in model)
			;        this is taken directly from module_mosaic_driver.F
			var_a = var_a * (180.0/211.0)


			;;;;; then extract the Oxygen data
			vbs_species = (/ "pcg1_f_o", "pcg2_f_o", "pcg3_f_o", "pcg4_f_o", \
							 "pcg5_f_o", "pcg6_f_o", "pcg7_f_o", "pcg8_f_o", "pcg9_f_o", \
							 "pcg1_b_o", "pcg2_b_o", "pcg3_b_o", "pcg4_b_o", \
							 "pcg5_b_o", "pcg6_b_o", "pcg7_b_o", "pcg8_b_o", "pcg9_b_o", \
							 "opcg1_f_o", "opcg2_f_o", "opcg3_f_o", "opcg4_f_o", \
							 "opcg5_f_o", "opcg6_f_o", "opcg7_f_o", "opcg8_f_o", \
							 "opcg1_b_o", "opcg2_b_o", "opcg3_b_o", "opcg4_b_o", \
							 "opcg5_b_o", "opcg6_b_o", "opcg7_b_o", "opcg8_b_o" \
							/)
			; obtain the oxygen mass
			var_b = vbs_species_extract_loop(a,species,vbs_string,vbs_species,inc_cw,time_in)

			;;;;; now calculate the O:C ratio
			; convert from ug/kg_air to ~mol/kg_air
			var_a = var_a / m_carbon
			var_b = var_b / m_oxygen
			
			var_a = var_b/where(var_a .gt. 1e-12, var_a, 1e-12)


			units = ""

		elseif (isStrSubset(species,"frac")) ; calculate mass fractions

			vbs_string = species
			
			;;;;; first extract the anthropogenic data
			vbs_species = (/ "pcg1_f_c", "pcg2_f_c", "pcg3_f_c", "pcg4_f_c", \
							 "pcg5_f_c", "pcg6_f_c", "pcg7_f_c", "pcg8_f_c", "pcg9_f_c", \
							 "pcg1_f_o", "pcg2_f_o", "pcg3_f_o", "pcg4_f_o", \
							 "pcg5_f_o", "pcg6_f_o", "pcg7_f_o", "pcg8_f_o", "pcg9_f_o", \
							 "opcg1_f_c", "opcg2_f_c", "opcg3_f_c", "opcg4_f_c", \
							 "opcg5_f_c", "opcg6_f_c", "opcg7_f_c", "opcg8_f_c", \
							 "opcg1_f_o", "opcg2_f_o", "opcg3_f_o", "opcg4_f_o", \
							 "opcg5_f_o", "opcg6_f_o", "opcg7_f_o", "opcg8_f_o" \
							/)



			vbs_species_2 = vbs_species
			;;; check which of these species are in the source file
			num_vbs = 0
			do ispec = 0,dimsizes(vbs_species)-1
				; check existence in source file
				if( isfilevar(a[0],vbs_species(ispec)+"_a01"))then
					vbs_species_2(num_vbs) = vbs_species(ispec)
					num_vbs = num_vbs + 1
				end if
			end do

			if(num_vbs .eq. 0)
				print("No anthropogenic VBS in file?!?")
			else
				; loop through the species listed for extraction
				spc_extract = vbs_species_2(0)+"_tot" ;str_sub_str(species, vbs_string, vbs_species_2(0))
				dummy = get_aerosol(a,spc_extract,inc_cw,time_in)
				var_a = (/ dummy /)
				desc = dummy@description
				delete(dummy)
				do ispec = 1, num_vbs-1
					spc_extract = vbs_species_2(ispec)+"_tot" ;str_sub_str(species, vbs_string, vbs_species_2(ispec))
					dummy = (/ get_aerosol(a,spc_extract,inc_cw,time_in) /)
					var_a = var_a + dummy
					delete(dummy)
				end do
			end if
			
			delete(vbs_species)
			delete(vbs_species_2)


			;;;;; then extract the biomass burning data
			vbs_species = (/ "pcg1_b_c", "pcg2_b_c", "pcg3_b_c", "pcg4_b_c", \
							 "pcg5_b_c", "pcg6_b_c", "pcg7_b_c", "pcg8_b_c", "pcg9_b_c", \
							 "pcg1_b_o", "pcg2_b_o", "pcg3_b_o", "pcg4_b_o", \
							 "pcg5_b_o", "pcg6_b_o", "pcg7_b_o", "pcg8_b_o", "pcg9_b_o", \
							 "opcg1_b_c", "opcg2_b_c", "opcg3_b_c", "opcg4_b_c", \
							 "opcg5_b_c", "opcg6_b_c", "opcg7_b_c", "opcg8_b_c", \
							 "opcg1_b_o", "opcg2_b_o", "opcg3_b_o", "opcg4_b_o", \
							 "opcg5_b_o", "opcg6_b_o", "opcg7_b_o", "opcg8_b_o" \
							/)

			vbs_species_2 = vbs_species
			;;; check which of these species are in the source file
			num_vbs = 0
			do ispec = 0,dimsizes(vbs_species)-1
				; check existence in source file
				if( isfilevar(a[0],vbs_species(ispec)+"_a01"))then
					vbs_species_2(num_vbs) = vbs_species(ispec)
					num_vbs = num_vbs + 1
				end if
			end do

			if(num_vbs .eq. 0)
				print("No biomass burning VBS in file?!?!")
				break
			else
				; loop through the species listed for extraction
				spc_extract = vbs_species_2(0)+"_tot" ;str_sub_str(species, vbs_string, vbs_species_2(0))
				dummy = get_aerosol(a,spc_extract,inc_cw,time_in)
				var_b = (/ dummy /)
				desc = dummy@description
				delete(dummy)
				do ispec = 1, num_vbs-1
					spc_extract = vbs_species_2(ispec)+"_tot" ;str_sub_str(species, vbs_string, vbs_species_2(ispec))
					dummy = (/ get_aerosol(a,spc_extract,inc_cw,time_in) /)
					var_b = var_b + dummy
					delete(dummy)
				end do
			end if

			delete(vbs_species)
			delete(vbs_species_2)

			

			;;;;; then extract the involatile core data
			vbs_species = (/ "oc" /)

			vbs_species_2 = vbs_species
			;;; check which of these species are in the source file
			num_vbs = 0
			do ispec = 0,dimsizes(vbs_species)-1
				; check existence in source file
				if( isfilevar(a[0],vbs_species(ispec)+"_a01"))then
					vbs_species_2(num_vbs) = vbs_species(ispec)
					num_vbs = num_vbs + 1
				end if
			end do

			if(num_vbs .eq. 0)
				print("No involatile OC in file?!?")
			else
				; loop through the species listed for extraction
				spc_extract = vbs_species_2(0)+"_tot" ;str_sub_str(species, vbs_string, vbs_species_2(0))
				dummy = get_aerosol(a,spc_extract,inc_cw,time_in)
				var_c = (/ dummy /)
				desc = dummy@description
				delete(dummy)
			end if

			
			;;; get total OC
			var_tot = var_a + var_b + var_c

			;;; determine if we're getting the fraction that is anthropogenic or biogenic
			if(isStrSubset(species,"anthro"))
			
				var_a = var_a/where(var_tot .gt. 1e-12, var_tot, 1e-12)

			elseif(isStrSubset(species,"bburn"))

				var_a = var_b/where(var_tot .gt. 1e-12, var_tot, 1e-12)

			elseif(isStrSubset(species,"invol"))

				var_a = var_c/where(var_tot .gt. 1e-12, var_tot, 1e-12)

			end if

			units = ""
			
			delete(var_tot)
			delete(var_b)
			delete(var_c)

		else ; not O:C solving or ratio calculation

			if(isStrSubset(species,"hoa"))
				vbs_string = "hoa"
				vbs_species = (/ "pcg1_f_c", "pcg2_f_c", "pcg3_f_c", "pcg4_f_c", \
								 "pcg5_f_c", "pcg6_f_c", "pcg7_f_c", "pcg8_f_c", "pcg9_f_c", \
								 "pcg1_f_o", "pcg2_f_o", "pcg3_f_o", "pcg4_f_o", \
								 "pcg5_f_o", "pcg6_f_o", "pcg7_f_o", "pcg8_f_o", "pcg9_f_o" \
								/)
			elseif(isStrSubset(species,"bboa"))
				vbs_string = "bboa"
				vbs_species = (/ "pcg1_b_c", "pcg2_b_c", "pcg3_b_c", "pcg4_b_c", \
								 "pcg5_b_c", "pcg6_b_c", "pcg7_b_c", "pcg8_b_c", "pcg9_b_c", \
								 "pcg1_b_o", "pcg2_b_o", "pcg3_b_o", "pcg4_b_o", \
								 "pcg5_b_o", "pcg6_b_o", "pcg7_b_o", "pcg8_b_o", "pcg9_b_o" \
								/)
			elseif(isStrSubset(species,"hsoa"))
				vbs_string = "hsoa"
				vbs_species = (/ "opcg1_f_c", "opcg2_f_c", "opcg3_f_c", "opcg4_f_c", \
								 "opcg5_f_c", "opcg6_f_c", "opcg7_f_c", "opcg8_f_c", \
								 "opcg1_f_o", "opcg2_f_o", "opcg3_f_o", "opcg4_f_o", \
								 "opcg5_f_o", "opcg6_f_o", "opcg7_f_o", "opcg8_f_o" \
								/)
			elseif(isStrSubset(species,"bbsoa"))
				vbs_string = "bbsoa"
				vbs_species = (/ "opcg1_b_c", "opcg2_b_c", "opcg3_b_c", "opcg4_b_c", \
								 "opcg5_b_c", "opcg6_b_c", "opcg7_b_c", "opcg8_b_c", \
								 "opcg1_b_o", "opcg2_b_o", "opcg3_b_o", "opcg4_b_o", \
								 "opcg5_b_o", "opcg6_b_o", "opcg7_b_o", "opcg8_b_o" \
								/)
			elseif(isStrSubset(species,"htoa"))
				vbs_string = "htoa"
				vbs_species = (/ "pcg1_f_c", "pcg2_f_c", "pcg3_f_c", "pcg4_f_c", \
								 "pcg5_f_c", "pcg6_f_c", "pcg7_f_c", "pcg8_f_c", "pcg9_f_c", \
								 "pcg1_f_o", "pcg2_f_o", "pcg3_f_o", "pcg4_f_o", \
								 "pcg5_f_o", "pcg6_f_o", "pcg7_f_o", "pcg8_f_o", "pcg9_f_o", \
								 "opcg1_f_c", "opcg2_f_c", "opcg3_f_c", "opcg4_f_c", \
								 "opcg5_f_c", "opcg6_f_c", "opcg7_f_c", "opcg8_f_c", \
								 "opcg1_f_o", "opcg2_f_o", "opcg3_f_o", "opcg4_f_o", \
								 "opcg5_f_o", "opcg6_f_o", "opcg7_f_o", "opcg8_f_o" \
								/)
			elseif(isStrSubset(species,"bbtoa"))
				vbs_string = "bbtoa"
				vbs_species = (/ "pcg1_b_c", "pcg2_b_c", "pcg3_b_c", "pcg4_b_c", \
								 "pcg5_b_c", "pcg6_b_c", "pcg7_b_c", "pcg8_b_c", "pcg9_b_c", \
								 "pcg1_b_o", "pcg2_b_o", "pcg3_b_o", "pcg4_b_o", \
								 "pcg5_b_o", "pcg6_b_o", "pcg7_b_o", "pcg8_b_o", "pcg9_b_o", \
								 "opcg1_b_c", "opcg2_b_c", "opcg3_b_c", "opcg4_b_c", \
								 "opcg5_b_c", "opcg6_b_c", "opcg7_b_c", "opcg8_b_c", \
								 "opcg1_b_o", "opcg2_b_o", "opcg3_b_o", "opcg4_b_o", \
								 "opcg5_b_o", "opcg6_b_o", "opcg7_b_o", "opcg8_b_o" \
								/)
			elseif(isStrSubset(species,"antsoa"))
				vbs_string = "antsoa"
				vbs_species = (/ "ant1_c", "ant2_c", "ant3_c", "ant4_c" \
								 /)
			elseif(isStrSubset(species,"biogsoa"))
				vbs_string = "biogsoa"
				vbs_species = (/ "biog1_c", "biog2_c", "biog3_c", "biog4_c" \
								/)
			elseif(isStrSubset(species,"soa"))
				vbs_string = "soa"
				vbs_species = (/ "ant1_c", "ant2_c", "ant3_c", "ant4_c", \
							 	 "biog1_c", "biog2_c", "biog3_c", "biog4_c" \
							 	/)
			elseif(isStrSubset(species,"totoa"))
				vbs_string = "totoa"
				vbs_species = (/ "pcg1_f_c", "pcg2_f_c", "pcg3_f_c", "pcg4_f_c", \
								 "pcg5_f_c", "pcg6_f_c", "pcg7_f_c", "pcg8_f_c", "pcg9_f_c", \
								 "pcg1_f_o", "pcg2_f_o", "pcg3_f_o", "pcg4_f_o", \
								 "pcg5_f_o", "pcg6_f_o", "pcg7_f_o", "pcg8_f_o", "pcg9_f_o", \
								 "pcg1_b_c", "pcg2_b_c", "pcg3_b_c", "pcg4_b_c", \
								 "pcg5_b_c", "pcg6_b_c", "pcg7_b_c", "pcg8_b_c", "pcg9_b_c", \
								 "pcg1_b_o", "pcg2_b_o", "pcg3_b_o", "pcg4_b_o", \
								 "pcg5_b_o", "pcg6_b_o", "pcg7_b_o", "pcg8_b_o", "pcg9_b_o", \
								 "opcg1_f_c", "opcg2_f_c", "opcg3_f_c", "opcg4_f_c", \
								 "opcg5_f_c", "opcg6_f_c", "opcg7_f_c", "opcg8_f_c", \
								 "opcg1_f_o", "opcg2_f_o", "opcg3_f_o", "opcg4_f_o", \
								 "opcg5_f_o", "opcg6_f_o", "opcg7_f_o", "opcg8_f_o", \
								 "opcg1_b_c", "opcg2_b_c", "opcg3_b_c", "opcg4_b_c", \
								 "opcg5_b_c", "opcg6_b_c", "opcg7_b_c", "opcg8_b_c", \
								 "opcg1_b_o", "opcg2_b_o", "opcg3_b_o", "opcg4_b_o", \
								 "opcg5_b_o", "opcg6_b_o", "opcg7_b_o", "opcg8_b_o", \
								 "ant1_c", "ant2_c", "ant3_c", "ant4_c", \
							     "biog1_c", "biog2_c", "biog3_c", "biog4_c", \
								 "oc" \
								/)
			elseif(isStrSubset(species,"OM"))
				vbs_string = "OM"
				vbs_species = (/ "pcg1_f_c", "pcg2_f_c", "pcg3_f_c", "pcg4_f_c", \
								 "pcg5_f_c", "pcg6_f_c", "pcg7_f_c", "pcg8_f_c", "pcg9_f_c", \
								 "pcg1_f_o", "pcg2_f_o", "pcg3_f_o", "pcg4_f_o", \
								 "pcg5_f_o", "pcg6_f_o", "pcg7_f_o", "pcg8_f_o", "pcg9_f_o", \
								 "pcg1_b_c", "pcg2_b_c", "pcg3_b_c", "pcg4_b_c", \
								 "pcg5_b_c", "pcg6_b_c", "pcg7_b_c", "pcg8_b_c", "pcg9_b_c", \
								 "pcg1_b_o", "pcg2_b_o", "pcg3_b_o", "pcg4_b_o", \
								 "pcg5_b_o", "pcg6_b_o", "pcg7_b_o", "pcg8_b_o", "pcg9_b_o", \
								 "opcg1_f_c", "opcg2_f_c", "opcg3_f_c", "opcg4_f_c", \
								 "opcg5_f_c", "opcg6_f_c", "opcg7_f_c", "opcg8_f_c", \
								 "opcg1_f_o", "opcg2_f_o", "opcg3_f_o", "opcg4_f_o", \
								 "opcg5_f_o", "opcg6_f_o", "opcg7_f_o", "opcg8_f_o", \
								 "opcg1_b_c", "opcg2_b_c", "opcg3_b_c", "opcg4_b_c", \
								 "opcg5_b_c", "opcg6_b_c", "opcg7_b_c", "opcg8_b_c", \
								 "opcg1_b_o", "opcg2_b_o", "opcg3_b_o", "opcg4_b_o", \
								 "opcg5_b_o", "opcg6_b_o", "opcg7_b_o", "opcg8_b_o", \
								 "ant1_c", "ant2_c", "ant3_c", "ant4_c", \
							     "biog1_c", "biog2_c", "biog3_c", "biog4_c", \
							     "oc" \
								/)
			elseif(isStrSubset(species,"POA"))
				vbs_string = "POA"
				vbs_species = (/ "pcg1_f_c", "pcg2_f_c", "pcg3_f_c", "pcg4_f_c", \
								 "pcg5_f_c", "pcg6_f_c", "pcg7_f_c", "pcg8_f_c", "pcg9_f_c", \
								 "pcg1_f_o", "pcg2_f_o", "pcg3_f_o", "pcg4_f_o", \
								 "pcg5_f_o", "pcg6_f_o", "pcg7_f_o", "pcg8_f_o", "pcg9_f_o", \
								 "pcg1_b_c", "pcg2_b_c", "pcg3_b_c", "pcg4_b_c", \
								 "pcg5_b_c", "pcg6_b_c", "pcg7_b_c", "pcg8_b_c", "pcg9_b_c", \
								 "pcg1_b_o", "pcg2_b_o", "pcg3_b_o", "pcg4_b_o", \
								 "pcg5_b_o", "pcg6_b_o", "pcg7_b_o", "pcg8_b_o", "pcg9_b_o", \
								 "opcg1_f_c", "opcg2_f_c", "opcg3_f_c", "opcg4_f_c", \
								 "opcg5_f_c", "opcg6_f_c", "opcg7_f_c", "opcg8_f_c", \
								 "opcg1_f_o", "opcg2_f_o", "opcg3_f_o", "opcg4_f_o", \
								 "opcg5_f_o", "opcg6_f_o", "opcg7_f_o", "opcg8_f_o", \
								 "opcg1_b_c", "opcg2_b_c", "opcg3_b_c", "opcg4_b_c", \
								 "opcg5_b_c", "opcg6_b_c", "opcg7_b_c", "opcg8_b_c", \
								 "opcg1_b_o", "opcg2_b_o", "opcg3_b_o", "opcg4_b_o", \
								 "opcg5_b_o", "opcg6_b_o", "opcg7_b_o", "opcg8_b_o", \
							     "oc" \
								/)
			elseif(isStrSubset(species,"vbs1"))
				vbs_string = "vbs1"
				vbs_species = (/ "pcg1_f_c", "pcg1_f_o", "opcg1_f_c", "opcg1_f_o", \
								 "pcg1_b_c", "pcg1_b_o", "opcg1_b_c", "opcg1_b_o" /)
			elseif(isStrSubset(species,"vbs2"))
				vbs_string = "vbs2"
				vbs_species = (/ "pcg2_f_c", "pcg2_f_o", "opcg2_f_c", "opcg2_f_o", \
								 "pcg2_b_c", "pcg2_b_o", "opcg2_b_c", "opcg2_b_o" /)
			elseif(isStrSubset(species,"vbs3"))
				vbs_string = "vbs3"
				vbs_species = (/ "pcg3_f_c", "pcg3_f_o", "opcg3_f_c", "opcg3_f_o", \
								 "pcg3_b_c", "pcg3_b_o", "opcg3_b_c", "opcg3_b_o" /)
			elseif(isStrSubset(species,"vbs4"))
				vbs_string = "vbs4"
				vbs_species = (/ "pcg4_f_c", "pcg4_f_o", "opcg4_f_c", "opcg4_f_o", \
								 "pcg4_b_c", "pcg4_b_o", "opcg4_b_c", "opcg4_b_o" /)
			elseif(isStrSubset(species,"vbs5"))
				vbs_string = "vbs5"
				vbs_species = (/ "pcg5_f_c", "pcg5_f_o", "opcg5_f_c", "opcg5_f_o", \
								 "pcg5_b_c", "pcg5_b_o", "opcg5_b_c", "opcg5_b_o" /)
			elseif(isStrSubset(species,"vbs6"))
				vbs_string = "vbs6"
				vbs_species = (/ "pcg6_f_c", "pcg6_f_o", "opcg6_f_c", "opcg6_f_o", \
								 "pcg6_b_c", "pcg6_b_o", "opcg6_b_c", "opcg6_b_o" /)
			elseif(isStrSubset(species,"vbs7"))
				vbs_string = "vbs7"
				vbs_species = (/ "pcg7_f_c", "pcg7_f_o", "opcg7_f_c", "opcg7_f_o", \
								 "pcg7_b_c", "pcg7_b_o", "opcg7_b_c", "opcg7_b_o" /)
			elseif(isStrSubset(species,"vbs8"))
				vbs_string = "vbs8"
				vbs_species = (/ "pcg8_f_c", "pcg8_f_o", "opcg8_f_c", "opcg8_f_o", \
								 "pcg8_b_c", "pcg8_b_o", "opcg8_b_c", "opcg8_b_o" /)
			elseif(isStrSubset(species,"vbs9"))
				vbs_string = "vbs9"
				vbs_species = (/ "pcg9_f_c", "pcg9_f_o", "pcg9_b_c", "pcg9_b_o" /)
			elseif(isStrSubset(species,"INORG"))
				vbs_string = "INORG"
				vbs_species = (/ "so4", "no3", "nh4", "na", "cl" \
								/)
			elseif(isStrSubset(species,"ALLDRY"))
				vbs_string = "ALLDRY"
				vbs_species = (/ "pcg1_f_c", "pcg2_f_c", "pcg3_f_c", "pcg4_f_c", \
								 "pcg5_f_c", "pcg6_f_c", "pcg7_f_c", "pcg8_f_c", "pcg9_f_c", \
								 "pcg1_f_o", "pcg2_f_o", "pcg3_f_o", "pcg4_f_o", \
								 "pcg5_f_o", "pcg6_f_o", "pcg7_f_o", "pcg8_f_o", "pcg9_f_o", \
								 "pcg1_b_c", "pcg2_b_c", "pcg3_b_c", "pcg4_b_c", \
								 "pcg5_b_c", "pcg6_b_c", "pcg7_b_c", "pcg8_b_c", "pcg9_b_c", \
								 "pcg1_b_o", "pcg2_b_o", "pcg3_b_o", "pcg4_b_o", \
								 "pcg5_b_o", "pcg6_b_o", "pcg7_b_o", "pcg8_b_o", "pcg9_b_o", \
								 "opcg1_f_c", "opcg2_f_c", "opcg3_f_c", "opcg4_f_c", \
								 "opcg5_f_c", "opcg6_f_c", "opcg7_f_c", "opcg8_f_c", \
								 "opcg1_f_o", "opcg2_f_o", "opcg3_f_o", "opcg4_f_o", \
								 "opcg5_f_o", "opcg6_f_o", "opcg7_f_o", "opcg8_f_o", \
								 "opcg1_b_c", "opcg2_b_c", "opcg3_b_c", "opcg4_b_c", \
								 "opcg5_b_c", "opcg6_b_c", "opcg7_b_c", "opcg8_b_c", \
								 "opcg1_b_o", "opcg2_b_o", "opcg3_b_o", "opcg4_b_o", \
								 "opcg5_b_o", "opcg6_b_o", "opcg7_b_o", "opcg8_b_o", \
								 "ant1_c", "ant2_c", "ant3_c", "ant4_c", \
								 "biog1_c", "biog2_c", "biog3_c", "biog4_c", \
								 "oc", "bc", "so4", "no3", "nh4", "na", "cl", "oin" \
								/)
			elseif(isStrSubset(species,"ALL"))
				vbs_string = "ALL"
				vbs_species = (/ "pcg1_f_c", "pcg2_f_c", "pcg3_f_c", "pcg4_f_c", \
								 "pcg5_f_c", "pcg6_f_c", "pcg7_f_c", "pcg8_f_c", "pcg9_f_c", \
								 "pcg1_f_o", "pcg2_f_o", "pcg3_f_o", "pcg4_f_o", \
								 "pcg5_f_o", "pcg6_f_o", "pcg7_f_o", "pcg8_f_o", "pcg9_f_o", \
								 "pcg1_b_c", "pcg2_b_c", "pcg3_b_c", "pcg4_b_c", \
								 "pcg5_b_c", "pcg6_b_c", "pcg7_b_c", "pcg8_b_c", "pcg9_b_c", \
								 "pcg1_b_o", "pcg2_b_o", "pcg3_b_o", "pcg4_b_o", \
								 "pcg5_b_o", "pcg6_b_o", "pcg7_b_o", "pcg8_b_o", "pcg9_b_o", \
								 "opcg1_f_c", "opcg2_f_c", "opcg3_f_c", "opcg4_f_c", \
								 "opcg5_f_c", "opcg6_f_c", "opcg7_f_c", "opcg8_f_c", \
								 "opcg1_f_o", "opcg2_f_o", "opcg3_f_o", "opcg4_f_o", \
								 "opcg5_f_o", "opcg6_f_o", "opcg7_f_o", "opcg8_f_o", \
								 "opcg1_b_c", "opcg2_b_c", "opcg3_b_c", "opcg4_b_c", \
								 "opcg5_b_c", "opcg6_b_c", "opcg7_b_c", "opcg8_b_c", \
								 "opcg1_b_o", "opcg2_b_o", "opcg3_b_o", "opcg4_b_o", \
								 "opcg5_b_o", "opcg6_b_o", "opcg7_b_o", "opcg8_b_o", \
								 "ant1_c", "ant2_c", "ant3_c", "ant4_c", \
							     "biog1_c", "biog2_c", "biog3_c", "biog4_c", \
								 "oc", "bc", "so4", "no3", "nh4", "na", "cl", "oin", "water" \
								/)
			elseif(isStrSubset(species,"DRY"))
				vbs_string = "DRY"
				vbs_species = (/ "pcg1_f_c", "pcg2_f_c", "pcg3_f_c", "pcg4_f_c", \
								 "pcg5_f_c", "pcg6_f_c", "pcg7_f_c", "pcg8_f_c", "pcg9_f_c", \
								 "pcg1_f_o", "pcg2_f_o", "pcg3_f_o", "pcg4_f_o", \
								 "pcg5_f_o", "pcg6_f_o", "pcg7_f_o", "pcg8_f_o", "pcg9_f_o", \
								 "pcg1_b_c", "pcg2_b_c", "pcg3_b_c", "pcg4_b_c", \
								 "pcg5_b_c", "pcg6_b_c", "pcg7_b_c", "pcg8_b_c", "pcg9_b_c", \
								 "pcg1_b_o", "pcg2_b_o", "pcg3_b_o", "pcg4_b_o", \
								 "pcg5_b_o", "pcg6_b_o", "pcg7_b_o", "pcg8_b_o", "pcg9_b_o", \
								 "opcg1_f_c", "opcg2_f_c", "opcg3_f_c", "opcg4_f_c", \
								 "opcg5_f_c", "opcg6_f_c", "opcg7_f_c", "opcg8_f_c", \
								 "opcg1_f_o", "opcg2_f_o", "opcg3_f_o", "opcg4_f_o", \
								 "opcg5_f_o", "opcg6_f_o", "opcg7_f_o", "opcg8_f_o", \
								 "opcg1_b_c", "opcg2_b_c", "opcg3_b_c", "opcg4_b_c", \
								 "opcg5_b_c", "opcg6_b_c", "opcg7_b_c", "opcg8_b_c", \
								 "opcg1_b_o", "opcg2_b_o", "opcg3_b_o", "opcg4_b_o", \
								 "opcg5_b_o", "opcg6_b_o", "opcg7_b_o", "opcg8_b_o", \
								 "ant1_c", "ant2_c", "ant3_c", "ant4_c", \
							     "biog1_c", "biog2_c", "biog3_c", "biog4_c", \
								 "oc", "bc", "so4", "no3", "nh4", "na", "cl", "oin" \
								/)
			end if

			vbs_species_2 = vbs_species
			;;; check which of these species are in the source file
			num_vbs = 0
			do ispec = 0,dimsizes(vbs_species)-1
				; check existence in source file
				if( isfilevar(a[0],vbs_species(ispec)+"_a01"))then
					vbs_species_2(num_vbs) = vbs_species(ispec)
					num_vbs = num_vbs + 1
				end if
			end do

			if(num_vbs .eq. 0)
				print("No VBS species in file!!!")
				break
			end if
			
			; record number of vbs species to be extracted
			;num_vbs = dimsizes(vbs_species)

			; loop through the species listed for extraction
			spc_extract = str_sub_str(species, vbs_string, vbs_species_2(0))
			dummy = get_aerosol(a,spc_extract,inc_cw,time_in)
			var_a = (/ dummy /)
			desc = dummy@description
			delete(dummy)
			do ispec = 1, num_vbs-1
				spc_extract = str_sub_str(species, vbs_string, vbs_species_2(ispec))
				dummy = (/ get_aerosol(a,spc_extract,inc_cw,time_in) /)
				var_a = var_a + dummy
				delete(dummy)
			end do
			
			; replace var name in description string, and add to data variable
			desc = str_sub_str(desc, vbs_species_2(0), vbs_string)
			var_a@description = desc
		
		end if


	else
          print("not aerosol species: needs '_' in name")
          break
	end if


     var_a@units = units
     return(var_a)
              
end

; FUNCTION: calc_ang_exp_loge
; Calculates the Angstrom Exponent between two known wavelengths,
; Given corresponding AODs.
; INPUTS:
;    * lambda_1 : first known wavelength (nm), must be single value
;    * tau_1    : AOD corresponding to wavelength 1.
;    * lambda_2 : Second known wavelength (nm), must be single value
;    * tau_2    : AOD corresponding to wavelength 2.

undef("calc_ang_exp_loge")
function calc_ang_exp_loge(lambda_1[1]:numeric,tau_1:numeric,\
                      lambda_2[1]:numeric,tau_2:numeric)
local tau_1_local, tau_2_local, lambda_1_local, lambda_2_local
begin
	tau_1_local = where(tau_1.gt.0.and.tau_2.gt.0,tau_1,1.0)
	tau_2_local = where(tau_1.gt.0.and.tau_2.gt.0,tau_2,1.0)
	lambda_1_local = where(lambda_1.gt.0.and.lambda_2.gt.0,lambda_1,1.0)
	lambda_2_local = where(lambda_1.gt.0.and.lambda_2.gt.0,lambda_2,1.0)
	
	Ang = ( log(tau_1_local/tau_2_local) / log(lambda_1_local/lambda_2_local) )
	return(Ang)
end



; FUNCTION: calc_new_tau_PNNL
; Estimates the AOD between two known AODs by calculating the angstrom exponent,
;   then applying it to another known AOD 
;
;   This follows the code in "module_aer_opt_out.F" - which, in turn, follows PNNL's lead
;
; INPUTS:
;    * lambda_1 : first known wavelength (nm), must be single value
;    * tau_1    : AOD corresponding to wavelength 1.
;    * lambda_2 : Second known wavelength (nm), must be single value
;    * tau_2    : AOD corresponding to wavelength 2.
;    * lambda_3 : Third known wavelength (nm), must be single value
;    * tau_3    : AOD corresponding to wavelength 3.
;    * lambda_new: wavelength corresponding to new AOD to be calculated.
;                  Must be between lambda_1 and lambda_3.

undef("calc_new_tau_PNNL")
function calc_new_tau_PNNL(lambda_1[1]:numeric,tau_1:numeric,\
                      		lambda_2[1]:numeric,tau_2:numeric,\
                      		lambda_3[1]:numeric,tau_3:numeric,\
                      		lambda_new[1]:numeric)
local ang
begin

	if ( ((lambda_new .gt. lambda_1) .and. (lambda_new .lt. lambda_3)) .or. \
		 ((lambda_new .lt. lambda_1) .and. (lambda_new .gt. lambda_3)) ) 	

		; first calculate angstrom coefficient
		ang  = calc_ang_exp_loge(lambda_1, tau_1, lambda_3, tau_3)
				
		; Now use angstrom exponent to calculate AOD at new wavelength
		tau_new = tau_2 * (lambda_2/lambda_new)^(-ang)
		
		return(tau_new)
	
	else
		print("calc_new_tau: new lambda must be between lambda_1 and lambda_3")
		break	
	end if
end



undef("calculate_aerosol_massfraction")
function calculate_aerosol_massfraction(FILES:list,var_name:string)	
local tail, tail_text, \
	var_temp_cl, var_temp_na, var_temp_no3, var_temp_so4, \
	var_temp_nh4, var_temp_oin, var_temp_bc, var_temp_om, \
	var_temp_total, \
	desc

begin

	; pull out the information about what aerosol we're processing
	if(isStrSubset(var_name,"_pm1"))
		tail = "_pm1"
		tail_text = "PM1 "
	elseif(isStrSubset(var_name,"_pm25"))
		tail = "_pm25"
		tail_text = "PM2.5 "
	elseif(isStrSubset(var_name,"_tot"))
		tail = "_tot"
		tail_text = "PM10 "
	else
		print("don't recognise this variable array '"+var_name+"', so not calculating aerosol mass fraction")
		break
	end if

		
	; load individual chemical species
	var_temp_cl  = get_aerosol(FILES,"cl"+tail,True,-1)	
	var_temp_na  = get_aerosol(FILES,"na"+tail,True,-1)	
	var_temp_no3 = get_aerosol(FILES,"no3"+tail,True,-1)	
	var_temp_so4 = get_aerosol(FILES,"so4"+tail,True,-1)	
	var_temp_nh4 = get_aerosol(FILES,"nh4"+tail,True,-1)	
	var_temp_oin = get_aerosol(FILES,"oin"+tail,True,-1)	
	var_temp_bc  = get_aerosol(FILES,"bc"+tail,True,-1)
	var_temp_om  = get_vbs_aerosol(FILES,"OM"+tail,True,-1)

	; combine total aerosol
	var_temp_total = var_temp_cl + var_temp_na + var_temp_no3 + \
					var_temp_so4 + var_temp_nh4 + var_temp_oin + \
					var_temp_om + var_temp_bc
	
	
	var_temp = new(dimsizes(var_temp_total),float,-9999)
	; determine which species we're getting the mass fraction for
	if(isStrSubset(var_name,"cl"))
		var_temp = where(var_temp_total.gt.0,var_temp_cl/var_temp_total,-9999)
		desc = tail_text+"Chloride mass fraction"
	elseif(isStrSubset(var_name,"na"))
		var_temp = where(var_temp_total.gt.0,var_temp_na/var_temp_total,-9999)
		desc = tail_text+"Sodium mass fraction"
	elseif(isStrSubset(var_name,"no3"))
		var_temp = where(var_temp_total.gt.0,var_temp_no3/var_temp_total,-9999)
		desc = tail_text+"Nitrate mass fraction"
	elseif(isStrSubset(var_name,"so4"))
		var_temp = where(var_temp_total.gt.0,var_temp_so4/var_temp_total,-9999)
		desc = tail_text+"Sulphate mass fraction"
	elseif(isStrSubset(var_name,"nh4"))
		var_temp = where(var_temp_total.gt.0,var_temp_nh4/var_temp_total,-9999)
		desc = tail_text+"Ammonium mass fraction"
	elseif(isStrSubset(var_name,"oin"))
		var_temp = where(var_temp_total.gt.0,var_temp_oin/var_temp_total,-9999)
		desc = tail_text+"Other inorganics mass fraction"
	elseif(isStrSubset(var_name,"bc"))
		var_temp = where(var_temp_total.gt.0,var_temp_bc/var_temp_total,-9999)
		desc = tail_text+"Black carbon mass fraction"
	elseif(isStrSubset(var_name,"OM"))
		var_temp = where(var_temp_total.gt.0,var_temp_om/var_temp_total,-9999)
		desc = tail_text+"OM mass fraction"
	else
		print("don't recognise this chemical variable '"+var_name+"', so not calculating aerosol mass fraction")
		break
	
	end if
	
	; Copy meta data from dummy	
	copy_VarMeta(var_temp_bc, var_temp) 

	var_temp@units       = ""
	var_temp@description = desc

	
	return(var_temp)

end



undef("calculate_aerosol_totalmass")
function calculate_aerosol_totalmass(FILES:list,var_name:string,vbs:logical)	
local tail, tail_text, \
	var_temp_cl, var_temp_na, var_temp_no3, var_temp_so4, \
	var_temp_nh4, var_temp_oin, var_temp_bc, var_temp_om, \
	var_temp

begin

	; pull out the information about what aerosol we're processing
	if(isStrSubset(var_name,"_pm1"))
		tail = "_pm1"
		tail_text = "PM1 total mass"
	elseif(isStrSubset(var_name,"_pm25"))
		tail = "_pm25"
		tail_text = "PM2.5 total mass"
	elseif(isStrSubset(var_name,"_tot"))
		tail = "_tot"
		tail_text = "PM10 total mass"
	else
		print("don't recognise this variable array '"+var_name+"', so not calculating aerosol mass fraction")
		break
	end if

		
	; load individual chemical species
	var_temp_cl  = get_aerosol(FILES,"cl"+tail,True,-1)	
	var_temp_na  = get_aerosol(FILES,"na"+tail,True,-1)	
	var_temp_no3 = get_aerosol(FILES,"no3"+tail,True,-1)	
	var_temp_so4 = get_aerosol(FILES,"so4"+tail,True,-1)	
	var_temp_nh4 = get_aerosol(FILES,"nh4"+tail,True,-1)	
	var_temp_oin = get_aerosol(FILES,"oin"+tail,True,-1)	
	var_temp_bc  = get_aerosol(FILES,"bc"+tail,True,-1)
	; check if we are processing a VBS aerosol species or not
	if(vbs)
		var_temp_om  = get_vbs_aerosol(FILES,"OM"+tail,True,-1)
	else
		var_temp_om  = get_aerosol(FILES,"oc"+tail,True,-1)
	end if


	; combine total aerosol
	var_temp = var_temp_cl + var_temp_na + var_temp_no3 + \
				var_temp_so4 + var_temp_nh4 + var_temp_oin + \
				var_temp_om + var_temp_bc
	
		
	var_temp@description = tail_text
	var_temp@units = "ug/kg-dryair"

	
	return(var_temp)

end



undef("column_data_calculation")
function column_data_calculation(FILES:list,var:string)
local var_name, PH, PHB, height, grid_thickness, pres, tk, ppm_conv, var_temp, dims, var_temp_a
begin
		; pull out the variable name
		var_name = str_sub_str(var, "_column","" )


		; checks for data type, and then extract variable data
		if(isStrSubset(var_name,"_a0") .or. isStrSubset(var_name,"_cw0") .or. \
			isStrSubset(var_name,"_A0") .or. isStrSubset(var_name,"_CW0") .or. \
			isStrSubset(var_name,"_fine") .or. isStrSubset(var_name,"_coarse") .or. \
			isStrSubset(var_name,"_pm1") .or. isStrSubset(var_name,"_pm25") .or. isStrSubset(var_name,"_tot") )
		
			print("column calculation not set up for aerosol species yet, so can't process "+var)
			exit

		else
			var_temp = wrf_user_getvar(FILES,var_name,-1)
			
			if(isatt(var_temp,"units"))
				if(.not.isStrSubset(var_temp@units,"ppm") .and. .not.isStrSubset(var_temp@units,"PPM"))
					print("column calculation requires input with units of ppmv, not "+var_temp@units)
					print("  you will need to modify the code in variable_loading_library to process "+var_name)
					exit
				end if
			else
				print("no units information for "+var_name+" - column calculation requires ppmv, ")
				print("   so we will assume this is what the units are, but please check to make sure this is right!")
			end if
			
		end if

		; pull out pressure height information, and calculate grid level thickness (in m)
		PH = wrf_user_getvar(FILES,"PH",-1)
		PHB = wrf_user_getvar(FILES,"PHB",-1)
		height = wrf_user_getvar(FILES,"height",-1)
		grid_thickness = get_grid_thickness(PH,PHB,height)


		; convert from ppm to molecules / cm3
		pres = wrf_user_getvar(FILES,"pres",-1)
		tk   = wrf_user_getvar(FILES,"tk",-1)
		ppm_conv  = ppm_moleconc_conv(pres, tk)
		var_temp = var_temp * ppm_conv

		; calculate molecules / cm2 / column for each grid cell (with conversion from m to cm for grid thickness)
		var_temp = var_temp * grid_thickness * 100.
		
		; sum molecules / cm2 for all columns
		dims = dimsizes(var_temp)
		if(dimsizes(dims).eq.4)
			var_temp_a = dim_sum_n_Wrap(var_temp,1)	
		elseif(dimsizes(dims).eq.3)
			var_temp_a = dim_sum_n_Wrap(var_temp,0)			
		end if
		var_temp_a@units = "molecules/cm2"

		return(var_temp_a)

end



;;;;;;; general variable loading routine ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

undef("load_variable")
function load_variable(FILES:list,var:string,vbs:logical,var_units:logical)
local var1, var2, var3, var_tail, var_temp, unit_old, unit_new, \
		lambda_1, tau_1, lambda_2, tau_2, lambda_3, tau_3, lambda_new, \
		dims, var_RH, var_cloud, tt, ii, jj, kk, \
		var_temp_a, var_temp_b, var_name, q_min, \
		w_a01, w_a02, w_a03, w_a04, w_a05, w_a06, w_a07, w_a08, \
		T, P, PB, tk, QVAPOR, QVAPOR_alt, rh_main, rh_reduced
begin

	
	; calculate NOx or NOy mixing ratios
	if(var .eq. "nox" .or. var .eq. "noy")
	
		no = wrf_user_getvar(FILES,"no",-1)
		no2 = wrf_user_getvar(FILES,"no2",-1)
	
		if(var .eq. "noy")

			hno3    = (/ wrf_user_getvar(FILES,"hno3",-1) /)
			no3     = (/ wrf_user_getvar(FILES,"no3",-1)  /)
			n2o5    = (/ wrf_user_getvar(FILES,"n2o5",-1) /)
			pan     = (/ wrf_user_getvar(FILES,"pan",-1)  /)
			hono    = (/ wrf_user_getvar(FILES,"hono",-1) /)
			hno4    = (/ wrf_user_getvar(FILES,"hno4",-1) /)
			onit    = (/ wrf_user_getvar(FILES,"onit",-1) /)
			nap     = (/ wrf_user_getvar(FILES,"nap",-1)  /)
			isopn   = (/ wrf_user_getvar(FILES,"isopn",-1)/)
        
			var_temp = (/ no + no2 + hno3 + no3 + 2.0*n2o5 + pan + \
							hono + hno4 + onit + nap + isopn /)
		
			copy_VarMeta (no, var_temp) 
			var_temp@description = "Mixing ratio of NOy"
 				
		else
			var_temp = (/no + no2/)
			copy_VarMeta (no, var_temp) 
			var_temp@description = "Mixing ratio of NOx"
		end if


	; If we want the Aerosol Optical Thickness at 550nm wavelength, then follow this process
	elseif(var .eq. "TAU55" .or. var .eq. "AOD550")
		;Calculate tau, from two given wavelengths (i.e. intermediate btwn the two lambda's)  
		lambda_1 = 300.                                ; wavelength1
		tau_1 = wrf_user_getvar(FILES,"TAUAER1",-1)         ; AOD data for wavelength1
		lambda_2 = 400.                                  ; wavelength2
		tau_2 = wrf_user_getvar(FILES,"TAUAER2",-1)         ; AOD data for wavelength2
		lambda_3 = 999.                                  ; wavelength2
		tau_3 = wrf_user_getvar(FILES,"TAUAER4",-1)         ; AOD data for wavelength2
		lambda_new = 550.                                ; new wavelength to be determined
	
		if(var .eq. "TAU55")
			var_temp = calc_new_tau_PNNL(lambda_1, tau_1, lambda_2, tau_2, lambda_3, tau_3, lambda_new)
		elseif(var .eq. "AOD550")
			var_temp_a = calc_new_tau_PNNL(lambda_1, tau_1, lambda_2, tau_2, lambda_3, tau_3, lambda_new)
			; sum over vertical axis (dimension #1)
			var_temp = dim_sum_n_Wrap(var_temp_a,1)
		end if

	; If we wish to get the column concentration of a gas through the atmosphere then follow this process
	elseif(isStrSubset(var,"_column"))

		
		var_temp = column_data_calculation(FILES,var)

	
	; If we wish to find the value of something below cloud then follow this process
	elseif(isStrSubset(var,"_belowcloud").or.isStrSubset(var,"_atcloud"))
	
		; define what our minimum QCLOUD (liquid water content) is for it to be a cloud
		; this is in g kg-1
		q_min = 0.01 ; defined low, to catch the edges of clouds
		; convert to kg kg-1
		q_min = q_min / 1000.
	
		if(isStrSubset(var,"_belowcloudbase"))
			var_name = str_sub_str(var, "_belowcloudbase","" )
		elseif(isStrSubset(var,"_belowcloud"))
			var_name = str_sub_str(var, "_belowcloud","" )
		elseif(isStrSubset(var,"_atcloudbase"))
			var_name = str_sub_str(var, "_atcloudbase","" )
		elseif(isStrSubset(var,"_atcloud"))
			var_name = str_sub_str(var, "_atcloud","" )
		else
			print("don't recognise variable ending: "+var)
			exit
		end if


		; extract variable data
		if(isStrSubset(var_name,"_a0") .or. isStrSubset(var_name,"_cw0") .or. \
			isStrSubset(var_name,"_A0") .or. isStrSubset(var_name,"_CW0") .or. \
			isStrSubset(var_name,"_fine") .or. isStrSubset(var_name,"_coarse") .or. \
			isStrSubset(var_name,"_pm1") .or. isStrSubset(var_name,"_pm25") .or. isStrSubset(var_name,"_tot") )
		
			; check if we are processing a VBS aerosol species or not
			if(.not.ismissing(vbs).and.vbs)
				var_temp_b = get_vbs_aerosol(FILES,var_name,True,-1)
			else
				var_temp_b = get_aerosol(FILES,var_name,True,-1)		
			end if

		elseif(var_name .eq. "depth")
			var_temp_b = wrf_user_getvar(FILES,"height",-1)
			dims = dimsizes(var_temp_b)
			var_temp_b(:,0:dims(1)-2,:,:) = var_temp_b(:,1:dims(1)-1,:,:) - var_temp_b(:,0:dims(1)-2,:,:)
			var_temp_b@description = "depth"
		
		elseif(var_name .eq. "RHreduceAW")
			;;; this is the calculation of the reduction in RH that would occur
			;;; if we removed aerosol water from water vapour.
		
			; load aerosol water variables
			w_a01 = wrf_user_getvar(FILES,"water_a01",-1)
			w_a02 = wrf_user_getvar(FILES,"water_a02",-1)
			w_a03 = wrf_user_getvar(FILES,"water_a03",-1)
			w_a04 = wrf_user_getvar(FILES,"water_a04",-1)
			w_a05 = wrf_user_getvar(FILES,"water_a05",-1)
			w_a06 = wrf_user_getvar(FILES,"water_a06",-1)
			w_a07 = wrf_user_getvar(FILES,"water_a07",-1)
			w_a08 = wrf_user_getvar(FILES,"water_a08",-1)
			; load atmospheric variables needed for RH
			T      = FILES[:]->T
			P      = FILES[:]->P
			PB     = FILES[:]->PB
			QVAPOR = FILES[:]->QVAPOR
			; process atmos variables (after WRFUserARW.ncl code for RH calcs)
			T = T + 300.
			P = P + PB
			tk = wrf_tk( P , T )
			; reduce QVAPOR (in kg/kg) by water in w_a0* (in ug/kg) - including unit conversion of 10^9
			QVAPOR_alt = QVAPOR - ( (w_a01+w_a02+w_a03+w_a04+w_a05+w_a06+w_a07+w_a08)/1e9 )
			; ensure we've no negative water vapour?!?!
			QVAPOR = QVAPOR > 0.000
			QVAPOR_alt = QVAPOR_alt > 0.000
			; calculate the new RH
			rh_main = wrf_rh( QVAPOR, P, tk )
			rh_reduced = wrf_rh( QVAPOR_alt, P, tk )
			
			; determine the reduction in RH from aerosol water
			var_temp_b = rh_main - rh_reduced
			
			var_temp_b@description = "RH reduction from aerosol water"
			var_temp_b@units = "%"
			
		else
			var_temp_b = wrf_user_getvar(FILES,var_name,-1)
		end if
		
		
		; extract cloud files, and check dimension sizes
		var_cloud = wrf_user_getvar(FILES,"QCLOUD",-1)
		dims = dimsizes(var_temp_b)
	
		; create array for storing output, and initialise as NaNs
		var_temp_a = new(dims,typeof(var_temp_b),-9999.)
		var_temp_a = -9999.
		
		
		if(isStrSubset(var,"belowcloud"))
		; Find the grid cells with no cloud, but with cloud in the grid cell above
		;
		; This will find all clear spaces below the bottom of clouds, not just at cloud base
		;
			var_temp_a(:,0:dims(1)-2,:,:) = \
						where( var_cloud(:,0:dims(1)-2,:,:).le.q_min .and. \
								var_cloud(:,1:dims(1)-1,:,:).gt.q_min , \
								var_temp_b(:,0:dims(1)-2,:,:), -9999)

			; if we are just interested in cloud base RH then we climb each column,
			;   looking for an RH value - recording it and exiting the loop once we find one 
			if(isStrSubset(var,"belowcloudbase"))

				var_temp = new((/dims(0),dims(2),dims(3)/),typeof(var_temp_b),-9999.)

				do tt = 0,dims(0)-1
				do ii = 0,dims(2)-1
				do jj = 0,dims(3)-1
					do kk = 0,dims(1)-1
						if(.not.ismissing(var_temp_a(tt,kk,ii,jj)))
							var_temp(tt,ii,jj) = var_temp_a(tt,kk,ii,jj)
							break
						end if
					end do
				end do
				end do
				end do
		
				var_temp@description = var_temp_b@description+" below cloud base"
		
			else
	
				var_temp = var_temp_a
				var_temp@description = var_temp_b@description+" below cloud"
	
			end if


		elseif(isStrSubset(var,"atcloudbase"))
		; Find the grid cells with cloud, but without cloud in the grid cell below
		;   (NOTE: this will ignore ground fog!)
		
			var_temp_a(:,1:dims(1)-1,:,:) = \
						where( var_cloud(:,0:dims(1)-2,:,:).le.q_min .and. \
								var_cloud(:,1:dims(1)-1,:,:).gt.q_min , \
								var_temp_b(:,1:dims(1)-1,:,:), -9999)
								
			var_temp = new((/dims(0),dims(2),dims(3)/),typeof(var_temp_b),-9999.)

			; now climb each column, selecting the first value, so giving the lowest 
			;  cloud base variable
			do tt = 0,dims(0)-1
			do ii = 0,dims(2)-1
			do jj = 0,dims(3)-1
				do kk = 0,dims(1)-1
					if(.not.ismissing(var_temp_a(tt,kk,ii,jj)))
						var_temp(tt,ii,jj) = var_temp_a(tt,kk,ii,jj)
						break
					end if
				end do
			end do
			end do
			end do
	
			var_temp@description = var_temp_b@description+" at cloud base"
			

		elseif(isStrSubset(var,"atcloud"))
		; Find the grid cells with cloud, but without cloud in the grid cell below

			var_temp = where(var_cloud.gt.q_min, var_temp_b, -9999)

			var_temp@description = var_temp_b@description+" in cloud"


		end if


		; copy the units over, if they exist
		if(isatt(var_temp_b,"units"))
			var_temp@units = var_temp_b@units
		end if
	
	
	; If we want to calculate the cloud water pH then we follow this process
	elseif(isStrSubset(var,"pH_")) 
	
		if (var .eq. "pH_cw")  
			var_temp = get_cloud_ph(FILES,-1)
		else ; assume this is an aerosol variable
			var_tail = str_sub_str(var, "pH_", "")
			var_temp = get_aerosol_ph(FILES,var_tail,-1)
		end if

	; If we want to calculate the wind speed
	elseif(isStrSubset(var,"windspeed")) 
	
		if (var .eq. "windspeed_10")  
			uu = wrf_user_getvar(FILES,"U10",-1)  ; u at 10 m
			vv = wrf_user_getvar(FILES,"V10",-1)  ; v at 10 m
		else ; assume we are getting the full volume windspeeds
			uu = wrf_user_getvar(FILES,"ua",-1)  
			vv = wrf_user_getvar(FILES,"va",-1)  
		end if
		
		var_temp = (uu^2.0 + vv^2.0)^0.5
		
	
	; If we want to calculate the wind direction
	elseif(isStrSubset(var,"winddirection")) 
	
		if (var .eq. "winddirection_10")  
			uvmet = wrf_user_getvar(FILES,"uvmet10",-1)  ; u & v at 10 m
			uu = uvmet(0,:,:,:)
			vv = uvmet(1,:,:,:)
		else ; assume we are getting the full volume winddirections
			uvmet = wrf_user_getvar(FILES,"uvmet",-1)  
			uu = uvmet(0,:,:,:,:)
			vv = uvmet(1,:,:,:,:)
		end if
		
		var_temp = new(dimsizes(uu),typeof(uu),-9999)
		
		; calculate wind direction
		var_temp = wind_direction(uu,vv,1)
		
		; set this to NaN where the wind is still
		var_temp = where(uu .ne. 0.0 .and. vv .ne. 0.0, var_temp, -9999)
	

	; calculate the aerosol mass fraction contribution from selected aerosol species
	elseif(isStrSubset(var,"_fraction"))
	
		; determine what species we're processing
		var_name_temp = str_sub_str(var, "_fraction", "")
		
		var_temp = calculate_aerosol_massfraction(FILES,var_name_temp)	
	
	
	; calculate the aerosol total mass
	elseif(isStrSubset(var,"totalmass_"))
	
		var_temp = calculate_aerosol_totalmass(FILES,var,vbs)	
	
	
	; Do quick and dirty check to see if input species name is an aerosol species.
	elseif(isStrSubset(var,"_a0") .or. isStrSubset(var,"_cw0") .or. \
			isStrSubset(var,"_A0") .or. isStrSubset(var,"_CW0") .or. \
			isStrSubset(var,"_fine") .or. isStrSubset(var,"_coarse") .or. \
			isStrSubset(var,"_pm1") .or. isStrSubset(var,"_pm25") .or. isStrSubset(var,"_tot") )
		
		; check if we are processing a VBS aerosol species or not
		if(.not.ismissing(vbs).and.vbs)
			var_temp = get_vbs_aerosol(FILES,var,True,-1)
		else
			var_temp = get_aerosol(FILES,var,True,-1)		
		end if
	
	else ; process as a gas species
		var_temp = wrf_user_getvar(FILES,var,-1)
	
	end if



	; check to see if we have specified a unit conversion for this variable
	if(isatt(var_units,var))
		unit_old = var_units@$var$(0)
		unit_new = var_units@$var$(1)
		print("converting from "+unit_old+" to "+unit_new)

		if(unit_old .eq. "ppmv") ;;;;;;;;;;;;;;; ppmv
			if(unit_new .eq. "ppbv")
				var_temp = var_temp * 1000.0
				var_temp@units = unit_new
			elseif(unit_new .eq. "pptv")
				var_temp = var_temp * 1000.0 * 1000.0
				var_temp@units = unit_new
			elseif(unit_new .eq. "ug/kg-dryair")
				if(var_units@$var$(2).eq."vbs")
					var_temp = var_temp * 1000.0  ;; first convert to ppbv
					var_temp = var_temp / (m_dry/m_vbs)
					var_temp@units = unit_new				
				else
					print("don't recognise molecular weight flag ("+var_units@$var$(2)+"), so not converting from ppmv to ug/kg-dryair")
				end if
			elseif(unit_new .eq. "molecules/cm3" .or. unit_new .eq. "molecules/cm~S~3~N~")
				; convert from ppm to molecules / cm3
				pres = wrf_user_getvar(FILES,"pres",-1)
				tk   = wrf_user_getvar(FILES,"tk",-1)
				ppm_conv  = ppm_moleconc_conv(pres, tk)
				var_temp = var_temp * ppm_conv
				var_temp@units = unit_new
			else
				print("don't recognise new unit '"+new_unit+"', so not converting units")		
			end if

		
		elseif(unit_old .eq. "K") ;;;;;;;;;;;;;; Kelvin
			if(unit_new .eq. "C")
				var_temp = var_temp - 273.15
				var_temp@units = unit_new
			else
				print("don't recognise new unit '"+new_unit+"', so not converting units")			
			end if

		elseif(unit_old .eq. "kg kg-1") ;;;;;;;;;;;;;; kg / kg (mixing ratio)
			if(unit_new .eq. "g kg-1")
				var_temp = var_temp * 1000.
				var_temp@units = unit_new
			elseif(unit_new .eq. "mg kg-1")
				var_temp = var_temp * 1000. * 1000.
				var_temp@units = unit_new
			elseif(unit_new .eq. "g scm-3")
				var_temp = var_temp * conv_kg_sm3 ; kg-1 to sm-3
				var_temp = var_temp * 1000. / 1e6  ; (kg->g) & (sm-3 -> scm-3)				
				var_temp@units = unit_new
			else
				print("don't recognise new unit '"+new_unit+"', so not converting units")			
			end if


		elseif(unit_old .eq. " kg-1") ;;;;;;;;;;;;;; # / kg (mixing ratio)
			if(unit_new .eq. "M# kg-1")
				var_temp = var_temp / ( 1000. * 1000. )
				var_temp@units = unit_new
			elseif(unit_new .eq. " scm-3")
				var_temp = var_temp * conv_kg_sm3 ; kg-1 to sm-3
				var_temp = var_temp / 1e6  ; sm-3 -> scm-3				
				var_temp@units = unit_new
			else
				print("don't recognise new unit '"+new_unit+"', so not converting units")			
			end if

		else
			print("don't recognise old unit '"+old_unit+"', so not converting units")
		
		end if



	end if


	return(var_temp)


end


;;;;;; file check routines
undef("check_file_status")
procedure check_file_status(Dfiles:string)
local itest, status, fail
begin

	;;; set fail flag to False
	fail = False

	;;; loop through the file list, making sure that they are all readable
	do itest = 0,dimsizes(Dfiles)-1
		status = systemfunc(" [ -r "+Dfiles(itest)+" ] && echo 'True' || echo 'False' ")
		;; if not readable then tell us, and set fail to True
		if(status .eq. "False")
			print(Dfiles(itest)+" is not readable")
			fail = True
		end if
	end do

	;; stop the program if the files aren't readable
	if(fail)
		print("Not all files are readable (see above messages) --- not continuing with script.")
		exit
	end if

end

;
;  Create a list of strings containing the paths to all the input files we wish to process.
;  This function assumes that the "domain" input lists only one domain - because processing 
;  more domains than this would break all our scripts. However the year, month, day, and hour
;  variables can contain multiple strings to process.
;
undef("create_input_file_list")
function create_input_file_list(plot_selected_times:logical,domain:string,\
				year_array:string,month_array:string,day_array:string,hour_array:string,\
				input_directory:string)
local File_List, ny, nm, nd, nh, iy, im, id, ih, name_string, Start_Array
begin
	

	; check if we need to filter the input files
	if(plot_selected_times)
		ny = dimsizes(year_array)
		nm = dimsizes(month_array)
		nd = dimsizes(day_array)
		nh = dimsizes(hour_array)
		
		if(ny.eq.1 .and. nm.eq.1 .and. nd.eq.1 .and. nh.eq.1)
			name_string = "wrfout_"+domain+"_"+year_array+"-"+month_array+"-"+day_array+"_"+hour_array+"*"
			File_List = systemfunc("ls -1 "+input_directory+"/"+name_string)	
		else
			Start_Array = True
			do iy = 0,ny-1
			do im = 0,nm-1
			do id = 0,nd-1
			do ih = 0,nh-1
				name_string = "wrfout_"+domain+"_"+year_array(iy)+"-"+month_array(im)+"-"+day_array(id)+"_"+hour_array(ih)+"*"
				print(" "+name_string)
				Temp_List = systemfunc("ls -1 "+input_directory+"/"+name_string)

				; check to see if we have only one file in the list
				if(dimsizes(Temp_List).eq.1)
					; check to make sure we have listed *some* files, and to jump to next iteration of 
					;   the loop if not
					if(ismissing(Temp_List))
						delete(Temp_List)
						continue
					end if
				end if
				
				
				
				if(Start_Array)
					File_List = Temp_List
					Start_Array = False
				else
					Old_List = File_List
					delete(File_List)
					File_List = array_append_record (Old_List, Temp_List, 0)
					delete(Old_List)
				end if
				delete(Temp_List)
			end do
			end do
			end do
			end do
			
			; check to see if we have failed to list any files - if stop program
			if(Start_Array)
				print("No input files for specified dates & times")
				exit
			end if
			
		end if
	else ; otherwise we just open everything for that domain
		File_List = systemfunc("ls -1 "+input_directory+"/wrfout_"+domain+"*")
	end if
	
	
	; sort the files into time order
	sqsort(File_List)
	
	print("Input Files:")
	do iy = 0,dimsizes(File_List)-1
		print(""+File_List(iy))
	end do

	return(File_List)

end 


;
;  creates input file list from array of file names
;
undef("create_multiple_input_file_list")
function create_multiple_input_file_list(plot_selected_times:logical,domain:string,\
				year_array:string,month_array:string,day_array:string,hour_array:string,\
				input_directory:string,number_scenarios:integer)
local File_List, il, Temp_List, Old_List, Start_Array
begin

	Start_Array = True

	do il = 0,number_scenarios-1
		Temp_List = create_input_file_list(plot_selected_times,"d01",year_array,month_array,day_array,hour_array,input_directory(il))

		if(Start_Array)
			File_List = Temp_List
			Start_Array = False
		else
			Old_List = File_List
			delete(File_List)
			File_List = array_append_record (Old_List, Temp_List, 0)
			delete(Old_List)
		end if
		delete(Temp_List)
		
	end do
	
	
	return(File_List)

end

;
;  Rough n' ready determination of the number of files in each scenario
;
;  This assumes that there's equal numbers of files in each. If there aren't
;   then you're doing something wrong!!!
;
undef("determine_start_end_positions")
function determine_start_end_positions(File_List:string,number_scenarios:integer)
local positions, num_files, num_files_scenario, ins
begin

	num_files = dimsizes(File_List)
	; check that there's no remainder
	if( (num_files % number_scenarios) .ne. 0 ) 
		print("there are not the same number of files in each scenario")
		print("total number of files: "+num_files)
		print("number of scenarios: "+number_scenarios)
		print("files per scenario: "+num_files/number_scenarios)
		exit 
	end if

	positions = new((/number_scenarios,2/),integer)

	num_files_scenario = num_files / number_scenarios

	do ins = 0,number_scenarios-1
		positions(ins,0) = num_files_scenario * ins
		positions(ins,1) = ( num_files_scenario * (ins+1) ) - 1
	end do

	return(positions)

end



;
;  Function for calculating the differences between two variables.
;  Note: The variables must have the exact same dimensions!
;
undef("difference_calculation")
function difference_calculation(plot_var_a:numeric,plot_var_b:numeric,percentage_difference:logical)
local plot_var
begin
	
	;;; absolute difference
	plot_var = where(.not.ismissing(plot_var_b).and. .not.ismissing(plot_var_a), plot_var_a - plot_var_b, -9999.0 )
	plot_var@_FillValue = -9999.0
	
	
	;;; calculate the percentage difference related to first input set (if requested)
	if(percentage_difference)
		plot_var = where(plot_var_a .gt. 0.0, plot_var / plot_var_a * 100.0, -9999.0) ; ensure we never divide by zero
		plot_var@units = "%"
	end if

	return(plot_var)

end


;
;  Function for creating sensible times strings to use for our file names.
;   This uses a function which has been created for setting up plot labels.
;
undef("create_sensible_time_strings")
function create_sensible_time_strings(FILES:list)
local time_strings, tempres, times, times_a, restick
begin

	;;; temporary resources variable to attach the attribute to
	tempres = True

	;;; get characters of the time information
	times = wrf_user_getvar(FILES,"Times",-1)
	;;; convert to integer string
	times_a = wrf_times_c(times,1)
	times_a@units = "hours since 1901-01-01 00:00:00"

	;;; create the time axis labels
	restick = True
	restick@ttmFormat = "%Y-%N-%D_%H"
	restick@ttmNumTicks = dimsizes(times_a)
	time_axis_labels(times_a,tempres,restick)	

	;;; copy the required strings to a new data array
	time_strings = tempres@tmXBLabels

	return(time_strings)

end


;
;  Function for creating times strings to use in our stats output files.
;   This uses a function which has been created for setting up plot labels.
;
undef("create_local_time_strings")
function create_local_time_strings(FILES:list,use_local_time:logical,istart:integer,iend:integer)
local time_strings, tempres, times_temp, times, times_a, restick
begin

	;;; temporary resources variable to attach the attribute to
	tempres = True

	;;; get characters of the time information
	times_temp = wrf_user_getvar(FILES,"Times",-1)
	times = times_temp(istart:iend,:)
	;;; convert to integer string
	times_a = wrf_times_c(times,1)
	times_a@units = "hours since 1901-01-01 00:00:00"

	;;; check if we are using local time, and if so change times (and label)
	if(use_local_time)
		times_a = times_a + use_local_time@offset
	end if

	;;; create the time axis labels
	restick = True
	restick@ttmFormat = "(%c-%D) %H"
	restick@ttmNumTicks = dimsizes(times_a)
	time_axis_labels(times_a,tempres,restick)	

	;;; copy the required strings to a new data array
	time_strings = tempres@tmXBLabels

	return(time_strings)

end





end  ;;; end of file
