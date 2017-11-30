;
; library file for the time series plotting scripts
;

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRFUserARW.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/contrib/time_axis_labels.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRF_contributed.ncl"

load "../../library/variable_loading_library.ncl"

begin

;;;;;;; Setting up the logical variables which we will use to carry 
;;;;;;;   attributes controlling the variable setups

;; used for carrying information about changing variable units, e.g.:
;;     variable_units@T2  = (/"K","C"/)
variable_units = True

;; used for carrying the variable descriptions
variable_description = True



;;; define if the species is a VBS compound or not
VBS_SPECIES = True
VBS_SPECIES@OM_pm1  = True
VBS_SPECIES@OM_pm25 = True
VBS_SPECIES@OM_tot  = True



;;;;;;; settings loading routine
undef("load_plot_settings")
procedure load_plot_settings(var:string,scenario:string,resources:logical,\
				use_local_time:logical) 
begin
	
	;;; delete the previous resources
	delete_VarAtts(resources,-1)
	
	;;; set general plotting resources
	resources@tmXBLabelAngleF = 90.
	if(use_local_time)
		if(use_local_time@offset .lt. 0)
			resources@tiXAxisString   = "Local Time (UTC "+use_local_time@offset+" )"
		else
			resources@tiXAxisString   = "Local Time (UTC +"+use_local_time@offset+" )"
		end if
	else
		resources@tiXAxisString   = "Time (UTC)"
	end if
	
end





undef("load_variable_and_sample_at_given_altitudes")
function load_variable_and_sample_at_given_altitudes(FILES:list,var:string,altitude_levels:numeric)
local var1, var2, var3, var_temp, variable, alt_km
begin

	;;; check if we are processing a VBS variable
	if(isatt(VBS_SPECIES,var))
		;;; call routine for extracting whole domain volume
		var_temp = load_variable(FILES,var,VBS_SPECIES@$var$,variable_units)		
	else
		;;; call routine for extracting whole domain volume, assume it's not a VBS species
		var_temp = load_variable(FILES,var,False,variable_units)
	end if

	;;; convert from metres to km for the altitude
	alt_km = altitude_levels / 1000.0


	;;; extract horizontal slices from the domain volume
	variable = wrf_user_vert_interp(FILES,var_temp,"ght_msl",alt_km,True)

	;;; replace any altitudes of zero with the data from the ground level
	do in = 0,dimsizes(alt_km)-1
		if(alt_km(in) .eq. 0.0)
			variable(:,in,:,:) = (/var_temp(:,0,:,:)/)
		end if
	end do



	return(variable)	
	
end


undef("load_2D_3D_variable_and_sample_at_given_altitudes")
function load_2D_3D_variable_and_sample_at_given_altitudes(FILES:list,var:string,altitude_levels:numeric)
local var1, var2, var3, var4, var5, var6, vara, dima, itime, var_temp, variable, alt_km
begin

	;;; check if we are processing a VBS variable
	if(isatt(VBS_SPECIES,var))
		;;; call routine for extracting whole domain volume
		var_temp = load_variable(FILES,var,VBS_SPECIES@$var$,variable_units)		

	;; process the 3D aerosol optical thickness data to give a total column aerosol optical depth
	else if(var .eq. "AOD55")
		var1 = load_variable(FILES,"TAU55",False,variable_units)
	
		; Find total column depth by summing through vertical levels:
		; declare variable to hold column data for variable. Set to ground level first
		var_temp = var1(:,0,:,:)	
		
		var_temp@_FillValue = -9999.
					
		dima = dimsizes(var1)
		
		do ik = 1, dima(1) - 1	; loop through vertical levels
			var_temp = var_temp + var1(:,ik,:,:)
		end do
		
	;; process the accumulated precipitation data to give hourly precip rates
	else if(isStrSubset(var,"_hour"))
		if(var .eq. "RAINC_hour")
			vara = load_variable(FILES,"RAINC",False,variable_units)	
		else if(var .eq. "RAINSH_hour")
			vara = load_variable(FILES,"RAINSH",False,variable_units)	
		else if(var .eq. "RAINNC_hour")
			vara = load_variable(FILES,"RAINNC",False,variable_units)	
		else if(var .eq. "SNOWNC_hour")
			vara = load_variable(FILES,"SNOWNC",False,variable_units)	
		else if(var .eq. "GRAUPELNC_hour")
			vara = load_variable(FILES,"GRAUPELNC",False,variable_units)	
		else if(var .eq. "HAILNC_hour")
			vara = load_variable(FILES,"HAILNC",False,variable_units)	
		else if(var .eq. "rain_hour" .or. var .eq. "precip_hour")
			var1 = load_variable(FILES,"RAINC",False,variable_units)	
			var2 = load_variable(FILES,"RAINSH",False,variable_units)	
			var3 = load_variable(FILES,"RAINNC",False,variable_units)	
			
			if(var .eq. "precip_hour")
				var4 = load_variable(FILES,"SNOWNC",False,variable_units)	
				var5 = load_variable(FILES,"GRAUPELNC",False,variable_units)	
				var6 = load_variable(FILES,"HAILNC",False,variable_units)	
				vara = var1+var2+var3+var4+var5+var6
			else	
				vara = var1+var2+var3
			end if
		else
			print("species "+var+" is unrecognised")
			break
		end if
		end if
		end if
		end if
		end if
		end if
		end if
	
		;; get the dimensions of the data set
		dima = dimsizes(vara)
	
		;; initialise the output variable
		var_temp = vara
		;; loop through the output variable, calculating the hourly precip rate
		do itime= 1,dima(0)-1
			var_temp(itime,:,:) = vara(itime,:,:) - vara(itime-1,:,:)
		end do
		;; set the first value of the output variable to zero
		var_temp(0,:,:) = 0.0 
		
		var_temp@units = "mm/hr"
	
		
	else ;;; extracting all other datasets
		;;; call routine for extracting whole domain volume, assume it's not a VBS species
		var_temp = load_variable(FILES,var,False,variable_units)
	end if
	end if
	end if


	;;;; sampling through the altitude levels *if* it's needed!

	;; first determine the dimensions of the model data
	dimensions = dimsizes(var_temp)
	; (time, height, x, y)  if data is 3D (plus time)
	;   or
	; (time, x, y)          if data is 2D (plus time)
	
	;; determine how many dimensions
	dimnum = dimsizes(dimensions)
	

	;;; if we have loaded a 3D variable then sample at given altitudes
	if ( dimnum .eq. 4 ) ; 3D data

		;;; convert from metres to km for the altitude
		alt_km = altitude_levels / 1000.0


		;;; extract horizontal slices from the domain volume
		variable = wrf_user_vert_interp(FILES,var_temp,"ght_msl",alt_km,True)

		;;; replace any altitudes of zero with the data from the ground level
		do in = 0,dimsizes(alt_km)-1
			if(alt_km(in) .eq. 0.0)
				variable(:,in,:,:) = (/var_temp(:,0,:,:)/)
			end if
		end do

	;;; otherwise just return the single level of data
	else if ( dimnum .eq. 3 ) ; 2D data

		variable = var_temp
		
	else
		print("don't know how to deal with data with "+dimsizes(dimensions)+" dimensions")
		exit
	end if
	end if
	


	return(variable)	
	
end




; calculates, standard stats for domain averages5th, 25th , median, 75th and 95th percentiles. 
; Returns them an array of size ( #time, #levels, 9):
;	stats(0) = mean(avg)
;	stats(1) = stdev
; 	stats(2) = min
;	stats(3) = max 
;	stats(4) = median
; 	stats(5) = lower quartile
;	stats(6) = upper quartile
;	stats(7) = 5th percentile
;	stats(8) = 95th percentile
undef("calculate_domain_stats")
function calculate_domain_stats(var_temp:numeric,domain_trim:logical)
local dimensions, dimnum, dimx, dimy, stat_data, \
		xstart, xend, ystart, yend, itime, ilev
begin
	;; first determine the dimensions of the model data
	dimensions = dimsizes(var_temp)
	; (time, height, x, y)  if data is 3D (plus time)
	;   or
	; (time, x, y)          if data is 2D (plus time)
	
	;; determine how many dimensions
	dimnum = dimsizes(dimensions)
	
	;; create stats storage array
	if ( dimnum .eq. 4 ) ; 3D data
		stat_data = new((/dimensions(0),dimensions(1),9/),float)
		dimx = dimensions(2)
		dimy = dimensions(3)
	else if ( dimnum .eq. 3 ) ; 2D data
		stat_data = new((/dimensions(0),9/),float)
		dimx = dimensions(1)
		dimy = dimensions(2)
	else
		print("don't know how to calculate stats for data with "+dimsizes(dimensions)+" dimensions")
		break
	end if
	end if
	
	;; check if we're trimming cells from the edge of the domain - and set the limits for the domain if we are
	if ( domain_trim )
		xstart = 0 + domain_trim@trim 
		xend   = dimx - 1 - domain_trim@trim
		ystart = 0 + domain_trim@trim 
		yend   = dimy - 1 - domain_trim@trim
	else
		xstart = 0 
		xend   = dimx - 1
		ystart = 0
		yend   = dimy - 1
	end if
	
	;; loop through the time dimension, calculating stats for each time slice
	do itime = 0,dimensions(0)-1
	
		if ( dimnum .eq. 4 )
			do ilev = 0,dimensions(1)-1
				stat_data(itime,ilev,:) = calc_standard_stats(var_temp(itime,ilev,xstart:xend,ystart:yend))
			end do
		else if ( dimnum .eq. 3 )
			stat_data(itime,:) = calc_standard_stats(var_temp(itime,xstart:xend,ystart:yend))
		end if
		end if
	
	end do
	
	;;; copy the units information
	if(isatt(var_temp,"units"))
		stat_data@units = var_temp@units
	end if

	return(stat_data)
	

end




;;;;;;; plotting routines ;;;;;;;;;;;;;;;;;;;;;;;


;---------------------------------------------------------------------;


; Original boxplot by  Adam Phillips
; Modified by Scott Archer-Nicholls
; Function for plotting box and whisker plots
; SAN 19/06/14: added new functionality:
;	Min and max of axis set by trYMinF and trYMaxF for vertical plots.
;	Min and max of axis set by trXMinF and trXMaxF for horizontal plots.

undef("boxplot_SAN")
function boxplot_SAN (wks:graphic,x[*]:numeric,y[*][*]:numeric,boxOpts:logical,plotres:logical,lineres:logical,HorizBox:logical)
;
;	This function creates a blank linLog plot object, on which box plots are created by extensive
;	use of gsn_add_polyline. The user can draw as many or as few box plots as they wish.
;	The plot is not drawn, and the frame is not advanced. May be used with gsn_panel.
;	The function options are as follows:
;
;	wks     ----- The workstation
;	x[*]    ----- A one-dimensional array containing the X-axis values of where the box plot(s)
;	              will be drawn.
;	y[*][*] ----- A two-dimensional array, where the rightmost dimension contains the box plot
;	 	      reference pt. data. y(n,0)=bottom_value, y(n,1)=bottom_value_of_box,
;		      y(n,2)=mid-value_of_box,y(n,3)=top_value_of_box,y(n,4)=top_value
;	boxOpts ----- 2 options attached as attributes may be attached here. 
;		      boxWidth ---- Sets the width of the box.
;		      boxColors ---- An array that contains the colors that the boxes will be drawn in.
;			             Ex. boxOpts@boxColors =  (/"green","blue"/)
;				     If the number of colors specified does not equal the number of 
;				     boxes drawn, only the first color specified will be used.
;	plotres ----- An optional xy-plot resource list. Will override defaults where applicable.
;	lineres ----- An optional resource list pertaining to the lines used to draw the boxes/lines.
;	!!! SAN - new option added:	!!!
;   HorizBox --- Logical, if True, will plot horizontally instead of vertically. It does
;                         this by reversing the x and y values, so the y-axis is defined 
;                         by x[*] and the x-axis by y[*][*].
;
;	Example call: plot3 = boxplot(wks,ntime,newyval,opti,res,False,False)
;		      draw(plot3)
;		      frame(wks)
;
local dimquery, numbox, min_axis_val, max_axis_val, xAxis, yAxis, dy, labarr,  boxwidth, \
	fsatts, yblab, polyres, boxcolor, cflag, gg, ff, xx, yy, dum, dumname, mres
begin
	dimquery = dimsizes(y)
	numbox = dimquery(0)
	if (numbox.ne.dimsizes(x)) then
	   print("Fatal: X must be one-dimensional and both X and Y must have the same rightmost dimension")
	   exit
	end if
	  

	if ( HorizBox) then

		; If we are doing horizontal plots, reverse all options so that x = y and y = x
		;	Developing y-axis
		yAxis = new(numbox+2,typeof(x))
		yAxis(1:numbox) = x
		if (numbox.ne.1) then
			dy = x(1)-x(0)
			yAxis(0) = x(0)-dy
			yAxis(numbox+1) = x(numbox-1)+dy
		else
			dy = 1
			yAxis(0) = x-dy
			yAxis(2) = x+dy
		end if
		
	
		if (boxOpts) then
			if (isatt(boxOpts,"boxWidth")) then
				boxwidth = boxOpts@boxWidth
			else
				boxwidth = dy*.3
			end if
		else
			boxwidth = dy*.3
		end if	
		
		labarr = new(numbox+2,"string")		;Prepare actual Y-axis labels...
		labarr(0) = ""
		labarr(numbox+1) = ""
		labarr(1:numbox) = yAxis(1:numbox)	

		; Whether to maximize plot in frame.
		maximize = get_res_value(plotres,"gsnMaximize",False) 
			
		; Check to see if the min and max on the axis has been set
		; If it is a horizontal plot, these should be set using the 
		; trXMinF and trXMaxF resources. 
		if (isatt(plotres,"trXMinF")) then
			min_axis_val= plotres@trXMinF
		else 
			min_axis_val =  min(y)-2
		end if

		if (isatt(plotres,"trXMaxF")) then
			max_axis_val = plotres@trXMaxF
		else 
			max_axis_val = max(y)+2
		end if 	


		if (plotres) then
			print("Plot resources detected, accepting")
			fsatts = getvaratts(plotres)
			do ty = 0,dimsizes(fsatts)-1
				if (fsatts(ty).eq."tmXBLabels") then		;Special section to test for XBLabels
					if (dimsizes(plotres@tmXBLabels).ne.numbox) then
						print("Fatal:Number of XB Labels does not match number of boxes, exiting")
						exit
					else
						labarr(1:numbox) = plotres@$fsatts(ty)$
					end if
					yblab = plotres@tmXBLabels
					delete(plotres@tmXBLabels)		;Delete so tmXBLabels is not used when all
				end if						;atts are assigned below...
			end do
			delete(fsatts)
		end if
		
		plot = create "plot" logLinPlotClass wks			
			"trXMinF" : min_axis_val
			"trXMaxF" : max_axis_val
			"trYMinF" : min(yAxis)
			"trYMaxF" : max(yAxis)
			"pmTickMarkDisplayMode" : "Always"
			"tmYLMode"   : "Explicit"
			"tmYLValues" : yAxis
			"tmYLLabels" : labarr
			"tmYROn"  : False
			"tmXTOn"  : False
			"tmYRBorderOn" : False
			"tmXTBorderOn" : False
			"pmTitleDisplayMode": "Always"                   ; allow titles
			"tiMainOn"     : True
			"tiMainString"   : ""	   
		end create
		if (plotres) then
			attsetvalues(plot,plotres)	
		end if

	; If HorizBox is false, do standard vertical box-plot
	else 

		;	Developing x-axis
		xAxis = new(numbox+2,typeof(x))
		xAxis(1:numbox) = x
		if (numbox.ne.1) then
			dx = x(1)-x(0)
			xAxis(0) = x(0)-dx
			xAxis(numbox+1) = x(numbox-1)+dx
		else
			dx = 1
			xAxis(0) = x-dx
			xAxis(2) = x+dx
		end if
	
		if (boxOpts) then
			if (isatt(boxOpts,"boxWidth")) then
				boxwidth = boxOpts@boxWidth
			else
				boxwidth = dx*.3
			end if
		else
			boxwidth = dx*.3
		end if	
		
		labarr = new(numbox+2,"string")		;Prepare actual X-axis labels...
		labarr(0) = ""
		labarr(numbox+1) = ""
		labarr(1:numbox) = xAxis(1:numbox)	

		; Whether to maximize plot in frame.
		maximize = get_res_value(plotres,"gsnMaximize",False) 
		
		; Check to see if the min and max on the axis has been set
		; If it is a vertical plot, these should be set using the 
		; trYMinF and trYMaxF resources. 
		if (isatt(plotres,"trYMinF")) then
			min_axis_val= plotres@trYMinF
		else 
			min_axis_val =  min(y)-2
		end if

		if (isatt(plotres,"trYMaxF")) then
			max_axis_val = plotres@trYMaxF
		else 
			max_axis_val = max(y)+2
		end if 	

		if (plotres) then
			print("Plot resources detected, accepting")
			fsatts = getvaratts(plotres)
			do ty = 0,dimsizes(fsatts)-1
				if (fsatts(ty).eq."tmXBLabels") then		;Special section to test for XBLabels
					if (dimsizes(plotres@tmXBLabels).ne.numbox) then
						print("Fatal:Number of XB Labels does not match number of boxes, exiting")
						exit
					else
						;labarr(1:numbox) = plotres@$fsatts(ty)$
						labarr(0:numbox-1) = plotres@$fsatts(ty)$
					end if
					xblab = plotres@tmXBLabels
					delete(plotres@tmXBLabels)		;Delete so tmXBLabels is not used when all
				end if					;atts are assigned below...
			end do
			delete(fsatts)
		end if
	
		plot = create "plot" logLinPlotClass wks
			"trYMinF" : min_axis_val
			"trYMaxF" : max_axis_val
			"trXMinF" : min(xAxis)
			"trXMaxF" : max(xAxis)
			"pmTickMarkDisplayMode" : "Always"
			"tmXBMode"   : "Explicit"
			"tmXBValues" : xAxis
			"tmXBLabels" : labarr
			"tmYROn"  : False
			"tmXTOn"  : False
			"tmYRBorderOn" : False
			"tmXTBorderOn" : False
			"pmTitleDisplayMode": "Always"                   ; allow titles
			"tiMainOn"     : True
			"tiMainString"   : ""	
		end create
		if (plotres) then
			attsetvalues(plot,plotres)	
		end if

	end if
	
	;--------------------------------------------------	
	; Section in middle here same in vertical and horizontal mode.
		
	polyres = True						;Set up defaults
	polyres@gsLineColor       = "black"                     ; color of lines
	polyres@gsLineThicknessF  = 1.5                       	; thickness of lines
	polyres@gsLineDashPattern = 0
	if (lineres) then
		fsatts = getvaratts(lineres)
		do ty = 0,dimsizes(fsatts)-1
			polyres@$fsatts(ty)$ = lineres@$fsatts(ty)$
		end do
	end if
	
	if (boxOpts) then
		if (isatt(boxOpts,"boxColors")) then
			boxcolor = boxOpts@boxColors
			if (dimsizes(boxcolor).eq.1.or.dimsizes(boxcolor).ne.numbox) then
				if (dimsizes(boxcolor).ne.numbox) then
					print("warning:Number of colors provided does not match number of boxes, using only one color")
				end if
				polyres@gsLineColor = boxcolor(0)
				cflag = 1
			else
				cflag = 2
			end if
		else
			cflag =1
		end if
	else
		cflag = 1
	end if
	
	;------------------------------------
	dum = new((/numbox,9/),graphic)
	
	; If in horizontal mode, reverse x and y.
	if ( HorizBox) then

		do gg = 0,numbox-1
		
		   ff = yAxis(gg+1)
		   
		   if (cflag.eq.2) then
			   polyres@gsLineColor = boxcolor(gg)
		   end if
		   
		   xx = (/y(gg,4),y(gg,4)/)
		   yy = (/(ff-(boxwidth/8.)),(ff+(boxwidth/8.))/)
		   dum(gg,0) = gsn_add_polyline(wks,plot,xx,yy,polyres)   

		   xx = (/y(gg,3),y(gg,4)/)
		   yy = (/ff,ff/)
		   polyres@gsLineDashPattern = 1
		   dum(gg,1) = gsn_add_polyline(wks,plot,xx,yy,polyres)    
		   polyres@gsLineDashPattern = 0

		   xx = (/y(gg,3),y(gg,3)/)
		   yy = (/(ff-(boxwidth/2.)),(ff+(boxwidth/2.))/)
		   dum(gg,2) = gsn_add_polyline(wks,plot,xx,yy,polyres)  

		   xx = (/y(gg,1),y(gg,3)/)
		   yy = (/(ff-(boxwidth/2.)),(ff-(boxwidth/2.))/)
		   dum(gg,3) = gsn_add_polyline(wks,plot,xx,yy,polyres)  

		   xx = (/y(gg,2),y(gg,2)/)
		   yy = (/(ff-(boxwidth/2.)),(ff+(boxwidth/2.))/)
		   dum(gg,4) = gsn_add_polyline(wks,plot,xx,yy,polyres)     

		   xx = (/y(gg,1),y(gg,3)/)
		   yy = (/(ff+(boxwidth/2.)),(ff+(boxwidth/2.))/)
		   dum(gg,5) = gsn_add_polyline(wks,plot,xx,yy,polyres)     

		   xx = (/y(gg,1),y(gg,1)/)
		   yy = (/(ff-(boxwidth/2.)),(ff+(boxwidth/2.))/)
		   dum(gg,6) = gsn_add_polyline(wks,plot,xx,yy,polyres)  

		   xx = (/y(gg,0),y(gg,1)/)
		   yy = (/ff,ff/)
		   polyres@gsLineDashPattern = 1
		   dum(gg,7) = gsn_add_polyline(wks,plot,xx,yy,polyres)  
		   polyres@gsLineDashPattern = 0

		   xx = (/y(gg,0),y(gg,0)/)
		   yy = (/(ff-(boxwidth/8.)),(ff+(boxwidth/8.))/)
		   dum(gg,8) = gsn_add_polyline(wks,plot,xx,yy,polyres)  
		end do
		dumname = unique_string("dum")								
		plot@$dumname$ = dum									
		if (isvar("yblab").eq."True") then							
		   plotres@tmXBLabels = yblab		;reassign XBLabels for possible future plots	
		end if				

	
	; If horiz is false, do standard vertical box-plot
	else 
	
		do gg = 0,numbox-1
		   ff = xAxis(gg+1)
		   if (cflag.eq.2) then
			   polyres@gsLineColor = boxcolor(gg)
		   end if
		   yy = (/y(gg,4),y(gg,4)/)
		   xx = (/(ff-(boxwidth/8.)),(ff+(boxwidth/8.))/)
		   dum(gg,0) = gsn_add_polyline(wks,plot,xx,yy,polyres)   

		   yy = (/y(gg,3),y(gg,4)/)
		   xx = (/ff,ff/)
		   polyres@gsLineDashPattern = 1
		   dum(gg,1) = gsn_add_polyline(wks,plot,xx,yy,polyres)   
		   polyres@gsLineDashPattern = 0

		   yy = (/y(gg,3),y(gg,3)/)
		   xx = (/(ff-(boxwidth/2.)),(ff+(boxwidth/2.))/)
		   dum(gg,2) = gsn_add_polyline(wks,plot,xx,yy,polyres)

		   yy = (/y(gg,1),y(gg,3)/)
		   xx = (/(ff-(boxwidth/2.)),(ff-(boxwidth/2.))/)
		   dum(gg,3) = gsn_add_polyline(wks,plot,xx,yy,polyres)   

		   yy = (/y(gg,2),y(gg,2)/)
		   xx = (/(ff-(boxwidth/2.)),(ff+(boxwidth/2.))/)
		   dum(gg,4) = gsn_add_polyline(wks,plot,xx,yy,polyres)   

		   yy = (/y(gg,1),y(gg,3)/)
		   xx = (/(ff+(boxwidth/2.)),(ff+(boxwidth/2.))/)
		   dum(gg,5) = gsn_add_polyline(wks,plot,xx,yy,polyres)   

		   yy = (/y(gg,1),y(gg,1)/)
		   xx = (/(ff-(boxwidth/2.)),(ff+(boxwidth/2.))/)
		   dum(gg,6) = gsn_add_polyline(wks,plot,xx,yy,polyres)   

		   yy = (/y(gg,0),y(gg,1)/)
		   xx = (/ff,ff/)
		   polyres@gsLineDashPattern = 1
		   dum(gg,7) = gsn_add_polyline(wks,plot,xx,yy,polyres)   
		   polyres@gsLineDashPattern = 0

		   yy = (/y(gg,0),y(gg,0)/)
		   xx = (/(ff-(boxwidth/8.)),(ff+(boxwidth/8.))/)
		   dum(gg,8) = gsn_add_polyline(wks,plot,xx,yy,polyres)
		end do
		dumname = unique_string("dum")								
		plot@$dumname$ = dum									
		if (isvar("xblab").eq."True") then							
		   plotres@tmXBLabels = xblab		;reassign XBLabels for possible future plots	
		end if				
	
	end if
	
	if(maximize) then
	   mres = True
	   mres@gsnDraw  = False
	   mres@gsnFrame = False
	   maximize_output(wks,mres)
    end if
	return(plot)
	
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;--------------------------------------------------------------------------;




;;;; stats plotting ;;;;;
undef("plot_stats_data")
procedure plot_stats_data(FILES:list,plot_stats:numeric,\
		var_name:string,output_directory:string,scenario:string,\
		var_resources:logical,\
		plot_type:string,altitude_levels:numeric,yaxis_limits:logical,\
		use_local_time:logical)
local plotres, times, times_a, restick, ntime, vardims, HorizBox, var_text, \
		output_file_name, output_file, wks, dimensions, dimscount, \
		newyval, ilev, plot3, plotres, lineres
	;;
	;;  plot_stats will enter with dimensions (time,[height],stats)   where stats has 9 values 
	;;
	;;   plotting script needs data in format:
	;	x[*]    ----- A one-dimensional array containing the X-axis values of where the box plot(s)
	;	              will be drawn.
	;	y[*][*] ----- A two-dimensional array, where the rightmost dimension contains the box plot
	;	 	      reference pt. data. y(n,0)=bottom_value, y(n,1)=bottom_value_of_box,
	;		      y(n,2)=mid-value_of_box,y(n,3)=top_value_of_box,y(n,4)=top_value
	;;
	;;
begin

	;;; copy plotting resources
	plotres = var_resources

	;;; get strings of the time information
	times = wrf_user_getvar(FILES,"Times",-1)
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
	time_axis_labels(times_a,plotres,restick)

	;;; set plotting locations for boxes, and set the x-axis labels and limits
	ntime = ispan(1,dimsizes(times_a),1)
	plotres@trXMinF  = min(times_a)-1
	plotres@trXMaxF  = max(times_a)+1

	;;; check if we want to set y-axis limits, and apply if we do
	if(yaxis_limits .and. isatt(yaxis_limits,var_name))
		plotres@trYMinF  = yaxis_limits@$var_name$(0)
		plotres@trYMaxF  = yaxis_limits@$var_name$(1)	
	end if

	;;; get dimensions of variable to plot
	vardims = dimsizes(plot_stats)
	
	;;; decide if we plot horizontally
	HorizBox = False

	;;; create the y-axis label
	if(isatt(variable_description,var_name))
		var_text = variable_description@$var_name$
	else
		var_text = var_name
	end if
	if(isatt(plot_stats,"units"))
		plotres@tiYAxisString = var_text + " (" +plot_stats@units+")"
	else
		plotres@tiYAxisString = var_text
	end if


	;;; name the output file
	output_file_name = var_name+"_"+scenario+"_domain_averages"
	output_file = output_directory+"/"+output_file_name
	

	;;; open plotting file
	wks = gsn_open_wks(plot_type,output_file)       ; Create a plot workstation



	;; determine how many dimensions to input data
	dimensions = dimsizes(plot_stats)
	dimscount = dimsizes(dimensions)
	
	;; create data array for the plotting data
	newyval = new((/dimensions(0),5/),float)


	;; determine box settings
	boxres = True
	boxres@boxColors = "blue"
	lineres = True
	lineres@gsLineThicknessF = 3.0


	if(dimscount .eq. 3) ; dealing with data that has altitude information
				
		do ilev = 0,dimensions(1)-1

			;;; create the plot title
			if(altitude_levels(ilev).eq.0.0)
				plotres@tiMainString = "Domain Statistics for "+var_text+" at ground level"
			else
				plotres@tiMainString = "Domain Statistics for "+var_text+" at "+flt2string(altitude_levels(ilev))+" metres ASL"
			end if
	
			newyval(:,0) = plot_stats(:,ilev,7)     ; 5th percentile
			newyval(:,1) = plot_stats(:,ilev,5)     ; lower quartile
			newyval(:,2) = plot_stats(:,ilev,4)     ; median
			newyval(:,3) = plot_stats(:,ilev,6)     ; upper quartile
			newyval(:,4) = plot_stats(:,ilev,8)     ; 95th percentile
	
		
			plot3 = boxplot_SAN(wks,times_a,newyval,boxres,plotres,lineres,HorizBox)
			draw(plot3)
			frame(wks)
	
	
		end do
	
	else if(dimscount .eq. 2) ; dealing with single level data

		;;; create the plot title
		plotres@tiMainString = "Domain Statistics for "+var_text
		
	
		newyval(:,0) = plot_stats(:,7)     ; 5th percentile
		newyval(:,1) = plot_stats(:,5)     ; lower quartile
		newyval(:,2) = plot_stats(:,4)     ; median
		newyval(:,3) = plot_stats(:,6)     ; upper quartile
		newyval(:,4) = plot_stats(:,8)     ; 95th percentile
	
		
		plot3 = boxplot_SAN(wks,times_a,newyval,boxres,plotres,lineres,HorizBox)
		draw(plot3)
		frame(wks)

	end if
	end if


	delete(wks)


end




;;;; stats plotting ;;;;;
undef("plot_single_location_timeseries")
procedure plot_single_location_timeseries(FILES:list,plot_var:numeric,\
		var_name:string,output_directory:string,scenario:string,\
		scenarios_array:string,var_resources:logical,\
		plot_type:string,plot_point:string,altitude_levels:numeric,\
		yaxis_limits:logical,use_local_time:logical,number_scenarios:integer,\
		scen_sep:integer)
local plotres, times_temp, times, times_a, times_a_plot, iscen, restick, ntime, plotres, loc, \
		vardims, var_text, output_file_name, output_file, wks, \
		dimensions, dimscount, ilev, plot3, var_scen
begin

	;;; copy plotting resources
	plotres = var_resources

	;;; get characters of the time information (and select just the first scenario information)
	times_temp = wrf_user_getvar(FILES,"Times",-1)
	times = times_temp(scen_sep(0,0):scen_sep(0,1),:)
	;;; convert to integer string
	times_a = wrf_times_c(times,1)
	times_a@units = "hours since 1901-01-01 00:00:00"
	
	;;; check if we are using local time, and if so change times (and label)
	if(use_local_time)
		times_a = times_a + use_local_time@offset
	end if

	;;; create time array for plotting
	times_a_plot = new((/number_scenarios,dimsizes(times_a)/),typeof(times_a))
	do iscen = 0,number_scenarios-1
		times_a_plot(iscen,:) = times_a(:)
	end do
	times_a_plot@units = times_a@units
	
	
	;;; create the time axis labels
	restick = True
	restick@ttmFormat = "(%c-%D) %H"
	restick@ttmNumTicks = dimsizes(times_a)
	time_axis_labels(times_a,plotres,restick)	
	
	;;; set the x-axis labels
	ntime = ispan(1,dimsizes(times_a),1)
	plotres@trXMinF  = min(times_a)
	plotres@trXMaxF  = max(times_a)

	;;; check if we want to set y-axis limits, and apply if we do
	if(yaxis_limits .and. isatt(yaxis_limits,var_name))
		plotres@trYMinF  = yaxis_limits@$var_name$(0)
		plotres@trYMaxF  = yaxis_limits@$var_name$(1)	
	end if


	;;; pull out the grid points for the plotted point
	loc = wrf_user_ll_to_ij(FILES[0], plot_point@Lon, plot_point@Lat, True)
	; convert to NCL subscripts
	loc = loc - 1


	;;; get dimensions of variable to plot
	vardims = dimsizes(plot_var)
	
	;;; create the y-axis label
	if(isatt(variable_description,var_name))
		var_text = variable_description@$var_name$
	else
		var_text = var_name
	end if
	if(isatt(plot_var,"units"))
		plotres@tiYAxisString = var_text + " (" +plot_var@units+")"
	else
		plotres@tiYAxisString = var_text
	end if


	;;; set the line plot resources
	plotres@xyMonoLineColor = False
	plotres@xyLineColors = (/"blue","green","red","cyan","purple"/)
	plotres@xyLineThicknessF = 3.0
	plotres@xyMarkLineMode = "MarkLines"
	plotres@xyMarkerColors = (/"blue","green","red","cyan","purple"/)
	plotres@xyMarkerSizes = (/10.0,10.0,10.0,10.0,10.0/)
	plotres@xyMarkers     = (/16,9,8,7,6/)



	;;; name the output file
	output_file_name = var_name+"_"+scenario+"_at_"+plot_point
	output_file = output_directory+"/"+output_file_name
	

	;;; open plotting file
	wks = gsn_open_wks(plot_type,output_file)       ; Create a plot workstation



	;; determine how many dimensions to input data
	dimensions = dimsizes(plot_var)
	dimscount = dimsizes(dimensions)
	

	if(dimscount .eq. 4) ; dealing with data that has altitude information
				
		do ilev = 0,dimensions(1)-1

			;;; create the plot title
			if(altitude_levels(ilev).eq.0.0)
				plotres@tiMainString = var_text+" at "+plot_point+" at ground level"
			else
				plotres@tiMainString = var_text+" at "+plot_point+" at "+flt2string(altitude_levels(ilev))+" metres ASL"
			end if
	
			;;; create the var data set for plotting
			var_scen = new((/number_scenarios,dimsizes(times_a)/),typeof(plot_var))
			do iscen = 0,number_scenarios-1
				var_scen(iscen,:) = plot_var(scen_sep(iscen,0):scen_sep(iscen,1),ilev,loc(1),loc(0))
			end do
	
			plot3 = gsn_csm_xy (wks,times_a_plot,var_scen,plotres)
				
			delete(var_scen)
	
		end do
		
	
	else if(dimscount .eq. 3) ; dealing with single level data

		;;; create the plot title
		plotres@tiMainString = var_name+" at "+plot_point

		;;; create the var data set for plotting
		var_scen = new((/number_scenarios,dimsizes(times_a)/),typeof(plot_var))
		do iscen = 0,number_scenarios-1
			var_scen(iscen,:) = plot_var(scen_sep(iscen,0):scen_sep(iscen,1),loc(1),loc(0))
		end do
		
		plot3 = gsn_csm_xy (wks,times_a_plot,var_scen,plotres)
			
		delete(var_scen)

	end if
	end if


	;***********************************************
	; legend resources
	;***********************************************
	labels = scenarios_array(0:number_scenarios-1)

	lgres                    = True

	lgres@vpWidthF           = 0.7                   ; width of legend (NDC)
	lgres@vpHeightF          = 0.4                   ; height of legend (NDC)

	lgres@lgTitleString      = ""
	lgres@lgTitleFontHeightF = 0.03
	lgres@lgAutoManage       = False                 ; Necessary to set font hgt

	lgres@lgMonoLineColor    = False
	lgres@lgMonoMarkerIndex  = False
	lgres@lgLineColors       = plotres@xyLineColors
	lgres@lgMonoItemType     = True                 ; more than one type
	lgres@lgItemType         = plotres@xyMarkLineMode   ; line/marker
	lgres@lgMarkerIndexes    = plotres@xyMarkers         ; marker indexes
	lgres@lgMarkerColors     = plotres@xyMarkerColors    ; marker colors
	lgres@lgMarkerSizes	     = plotres@xyMarkerSizes
	lgres@lgLineThicknessF   = plotres@xyLineThicknessF
	lgres@lgLineLabelStrings = (/"","","","",""/)        ; blank labels

	;
	; Draw the legend, indicating the number of items, a label for each
	; item, and the X, Y position of the legend in NDC coordinates.
	;
	gsn_legend_ndc(wks,number_scenarios,labels,0.1,0.95,lgres)

	frame(wks)                                        ; finally advance frame



end



;;; data writing scripts ;;;
undef("write_stats_data")
procedure write_stats_data(FILES:list,plot_stats:numeric,var_name:string,\
							output_directory:string,scenario:string,\
							altitude_levels:numeric,use_local_time:logical)
local dimensions, dimscount, time_strings, file_names, ilev, header_strings, \
		col_head_string, columns, itim, cf, stat_order, level_limit, out_string_array
begin

	;; determine how many dimensions to input data
	dimensions = dimsizes(plot_stats)
	dimscount = dimsizes(dimensions)
	
	; get the time strings for our data sets
	time_strings = create_local_time_strings(FILES,use_local_time,0,dimensions(0)-1)


	if(dimscount .eq. 3) ; dealing with data that has altitude information				

		; create the filenames and header strings
		file_names = new(dimsizes(altitude_levels),string)
		header_strings = new(dimsizes(altitude_levels),string)
		do ilev = 0,dimsizes(altitude_levels)-1
			if(altitude_levels(ilev).eq.0)
				file_names(ilev) = var_name+"_ground_level_domain_stats.txt"
				if(isatt(plot_stats,"units"))
					header_strings(ilev) = var_name+" ("+plot_stats@units+") at ground level"
				else
					header_strings(ilev) = var_name+" at ground level"				
				end if
			else
				file_names(ilev) = var_name+"_"+flt2string(altitude_levels(ilev))+"m_ASL_domain_stats.txt"
				if(isatt(plot_stats,"units"))
					header_strings(ilev) = var_name+" ("+plot_stats@units+") at "+flt2string(altitude_levels(ilev))+"m altitude (above sea level)"
				else
					header_strings(ilev) = var_name+" at "+flt2string(altitude_levels(ilev))+"m altitude (above sea level)"
				end if
			end if
		end do

		; set the number of levels for the rest of this subroutine
		level_limit = dimsizes(altitude_levels)


	else if(dimscount .eq. 2) ; dealing with single level data
	
		; create the filenames and header strings
		file_names = new(1,string)
		header_strings = new(1,string)

		file_names(0) = var_name+"_domain_stats.txt"
		if(isatt(plot_stats,"units"))
			header_strings(0) = var_name+" ("+plot_stats@units+")"
		else
			header_strings(0) = var_name
		end if

		; set the number of levels for the rest of this subroutine
		level_limit = 1

	end if
	end if
	

	; create the column header string
	if(use_local_time)
		col_head_string = "Times (UTC"+use_local_time@offset+"), Minimum, 5th Percentile, 25th Percentile, "+\
							"Median, 75th Percentile, 95th Percentile, Maximum, Mean"
	else
		col_head_string = "Times (UTC), Minimum, 5th Percentile, 25th Percentile, "+\
							"Median, 75th Percentile, 95th Percentile, Maximum, Mean"
	end if

	; create the column data
	columns = new((/dimsizes(time_strings),level_limit/),string)

	
	; set the delimiter for the columns
	cf = ", "
	; set the statistical data order
	stat_order = (/2, 7, 5, 4, 6, 8, 3, 0/)
	
	
	; loop through time and levels, to create strings
	do itim = 0, dimsizes(time_strings)-1
	do ilev = 0, level_limit-1
		columns(itim,ilev) = time_strings(itim)
		; loop through the stats, getting them into the order defined above
		do istat = 0,dimsizes(stat_order)-1
			if(dimscount .eq. 3) ; dealing with data that has altitude information				
				columns(itim,ilev) = columns(itim,ilev) + cf+sprintf("%7.4f",plot_stats(itim,ilev,stat_order(istat)))
			else if(dimscount .eq. 2) ; dealing with single level data
				columns(itim,ilev) = columns(itim,ilev) + cf+sprintf("%7.4f",plot_stats(itim,stat_order(istat)))
			end if
			end if
		end do
	end do
	end do

	; loop through the model levels, organise the data strings, and write them out
	do ilev = 0, level_limit-1
	    out_string_array = new(2+dimsizes(time_strings),string)
	    out_string_array(0) = header_strings(ilev)
	    out_string_array(1) = col_head_string
	    do istat = 2,dimsizes(out_string_array)-1	    
	    	out_string_array(istat) = columns(istat-2,ilev)
	    end do
	    asciiwrite(output_directory+"/"+file_names(ilev),out_string_array)
	    delete(out_string_array)
	end do



end




;;; data writing scripts ;;;
undef("write_single_location_timeseries")
procedure write_single_location_timeseries(FILES:list,plot_var:numeric,var_name:string,\
							output_directory:string,scenario:string,\
							scenarios_array:string,plot_point:string,\
							altitude_levels:numeric,use_local_time:logical,\
							number_scenarios:integer,scen_sep:integer)



local dimensions, dimscount, time_string_temp, time_strings, loc, loc_string, \
		file_names, header_strings, ilev, level_limit, \
		col_head_string, iscen, columns, cf, var_scen, \ 
		out_string_array, itim
begin
	
	;; determine how many dimensions to input data
	dimensions = dimsizes(plot_var)
	dimscount = dimsizes(dimensions)

	; get the time strings for our data sets (using just one set of data files for one scenario)
	time_strings = create_local_time_strings(FILES,use_local_time,scen_sep(0,0),scen_sep(0,1))

	;;; pull out the grid points for the plotted point
	loc = wrf_user_ll_to_ij(FILES[0], plot_point@Lon, plot_point@Lat, True)
	; convert to NCL subscripts
	loc = loc - 1

	;;; create location information string
	loc_string = plot_point+", Longitude "+flt2string(plot_point@Lon)+", Latitude "+flt2string(plot_point@Lat)


	if(dimscount .eq. 4) ; dealing with data that has altitude information				

		; create the filenames and header strings
		file_names = new(dimsizes(altitude_levels),string)
		header_strings = new(dimsizes(altitude_levels),string)
		do ilev = 0,dimsizes(altitude_levels)-1
			if(altitude_levels(ilev).eq.0)
				file_names(ilev) = var_name+"_at_"+plot_point+"_at_ground_level.txt"
				if(isatt(plot_var,"units"))
					header_strings(ilev) = var_name+" ("+plot_var@units+") at ground level, "+loc_string
				else
					header_strings(ilev) = var_name+" at ground level, "+loc_string
				end if
			else
				file_names(ilev) = var_name+"_at_"+plot_point+"_at_"+flt2string(altitude_levels(ilev))+"m_ASL.txt"
				if(isatt(plot_var,"units"))
					header_strings(ilev) = var_name+" ("+plot_var@units+") at "+flt2string(altitude_levels(ilev))+"m altitude (above sea level), "+loc_string
				else
					header_strings(ilev) = var_name+" at "+flt2string(altitude_levels(ilev))+"m altitude (above sea level), "+loc_string
				end if
			end if
		end do

		; set the number of levels for the rest of this subroutine
		level_limit = dimsizes(altitude_levels)


	else if(dimscount .eq. 3) ; dealing with single level data
	
		; create the filenames and header strings
		file_names = new(1,string)
		header_strings = new(1,string)

		file_names(0) = var_name+"_at_"+plot_point+".txt"
		if(isatt(plot_var,"units"))
			header_strings(0) = var_name+" ("+plot_var@units+"), "+loc_string
		else
			header_strings(0) = var_name+", "+loc_string
		end if

		; set the number of levels for the rest of this subroutine
		level_limit = 1

	end if
	end if
	

	; create the column header string
	if(use_local_time)
		col_head_string = "Times (UTC"+use_local_time@offset+")"
	else
		col_head_string = "Times (UTC)"
	end if
	do iscen = 0,number_scenarios-1
		col_head_string = col_head_string+", "+scenarios_array(iscen)
	end do

	; create the column data
	columns = new((/dimsizes(time_strings),level_limit/),string)
	

	; extract the var data sets for plotting
	var_scen = new((/number_scenarios,level_limit,dimsizes(time_strings)/),typeof(plot_var))
	if(dimscount .eq. 4) ; dealing with data that has altitude information				
		do ilev = 0, level_limit-1
		do iscen = 0, number_scenarios-1
				var_scen(iscen,ilev,:) = plot_var(scen_sep(iscen,0):scen_sep(iscen,1),ilev,loc(1),loc(0))
		end do
		end do
	else if(dimscount .eq. 3) ; dealing with single level data
		do iscen = 0,number_scenarios-1
				var_scen(iscen,0,:) = plot_var(scen_sep(iscen,0):scen_sep(iscen,1),loc(1),loc(0))
		end do
	end if
	end if


	; set the delimiter for the columns
	cf = ", "


	; loop through time and levels, to create strings which start with date information
	do itim = 0, dimsizes(time_strings)-1
	do ilev = 0, level_limit-1
		columns(itim,ilev) = time_strings(itim)
		; loop through the stats, getting them into the order defined above
		do iscen = 0,number_scenarios-1
			columns(itim,ilev) = columns(itim,ilev) + cf+sprintf("%7.4f",var_scen(iscen,ilev,itim))
		end do
	end do
	end do



	; loop through the model levels, organise the data strings, and write them out
	do ilev = 0, level_limit-1
	    out_string_array = new(2+dimsizes(time_strings),string)
	    out_string_array(0) = header_strings(ilev)
	    out_string_array(1) = col_head_string
	    do itim = 2,dimsizes(out_string_array)-1	    
	    	out_string_array(itim) = columns(itim-2,ilev)
	    end do
	    asciiwrite(output_directory+"/"+file_names(ilev),out_string_array)
	    delete(out_string_array)
	end do



end







end ;; end of script
