;
;  Statistical Functions  (from WRFCHEM_misc.ncl)
;



load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"     ; include standard resources
load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRFUserARW.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

begin





;--------------------------------------------------------------------------------

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; Functions for calculating statistics from data							;;;;;
;--------------------------------------------------------------------------------

; rm_missing_1d - remove missing values from a 1D array and resize array
undef("rm_missing_1d") 
function rm_missing_1d(data[*]:numeric)
local ii
begin
	ii = ind(.not.ismissing(data))
	oneD_data_processed = data(ii)
	return(oneD_data_processed)
end

;--------------------------------------------------------------------------------

; Calculates weighted mean.
; Requires two input arrays - must be of the same size (or will crash)
; Based on function:
; weighted_mean = sum(x_i*w_i)/sum(w_i)

undef("calc_wght_mean")
function calc_wght_mean(data:numeric, weight:numeric)
local Sigma_data, Sigma_wght, weighted_data
begin
	; WARNING: input array data must be of same size as Weight array

	Sigma_wght = sum(weight)

	weighted_data = data * weight
	Sigma_data = sum(weighted_data)

	wght_mean = (Sigma_data)/Sigma_wght

	return(wght_mean)

end

;--------------------------------------------------------------------------------

; Calculates median by sorting data into 1D array
undef("calc_median")
function calc_median(data:numeric)
local oneD_data 
begin
	oneD_data = ndtooned(data)		; put data into 1 dimension
	qsort(oneD_data)	; sort data - ascending order
	median = dim_median(oneD_data)

	return(median)
end

;--------------------------------------------------------------------------------
; Calculates lower quartile by sorting data into 1D array
undef("calc_lower_quartile")
function calc_lower_quartile(data:numeric)
local oneD_data, numel, pt_x25, oneD_data_in, ii
begin
	oneD_data_in = ndtooned(data)		; put data into 1 dimension

	ii = ind(.not.ismissing(oneD_data_in))  ; find locations of non-missing values
	oneD_data = oneD_data_in(ii)            ; new array containing no missing values


	qsort(oneD_data)	; sort data - ascending order
	numel = dimsizes(oneD_data)	; find number of elements
		
	; find point of 25th percentile
	pt_x25 = round(.25*(numel-1),3)
	lower_quart = oneD_data(pt_x25)
	
	return(lower_quart)
end
	
;--------------------------------------------------------------------------------
; similar for upper quartile
undef("calc_upper_quartile")
function calc_upper_quartile(data:numeric)
local oneD_data, sorted_data, numel, pt_x75, oneD_data_in, ii
begin
	oneD_data_in = ndtooned(data)		; put data into 1 dimension

	ii = ind(.not.ismissing(oneD_data_in))  ; find locations of non-missing values
	oneD_data = oneD_data_in(ii)            ; new array containing no missing values

	qsort(oneD_data)	; sort data - ascending order
	numel = dimsizes(oneD_data)	; find number of elements
		
	; find point of 75th percentile
	pt_x75 = round(.75*(numel-1),3)
	upper_quart = oneD_data(pt_x75)
	
	return(upper_quart)
end

;--------------------------------------------------------------------------------
; calculation for given percentile
undef("calc_percentile")
function calc_percentile(data:numeric, perc_point:numeric)
local oneD_data, sorted_data, numel, pt_x, oneD_data_in, ii
begin
	oneD_data_in = ndtooned(data)		; put data into 1 dimension

	ii = ind(.not.ismissing(oneD_data_in))  ; find locations of non-missing values
	oneD_data = oneD_data_in(ii)            ; new array containing no missing values

	qsort(oneD_data)				; sort data - ascending order
	numel = dimsizes(oneD_data)		; find number of elements
	
	; find point of percentile - perc_point must be between 0 and 100
	pt_x = round(.01*perc_point*(numel-1),3)
	percentile = oneD_data(pt_x)
	
	return(percentile)
end

;--------------------------------------------------------------------------------
; calc_boxplot_stats: combination of above functions
; calculates, 5th, 25th , median, 75th and 95th percentiles. 
; Returns them an array of size 5:
; 	stats(0) = 5th percentile
; 	stats(1) = lower quartile
;	stats(2) = median
;	stats(3) = upper quartile
;	stats(4) = 95th percentile

undef("calc_boxplot_stats")
function calc_boxplot_stats(data:numeric)
local oneD_data, numel, perc_point, pt_x, oneD_data_in, ii, j
begin

	ndims = dimsizes(dimsizes(data))
	if (ndims .gt. 1)
		oneD_data_in = ndtooned(data) 
	else
		oneD_data_in = data 
	end if
	
	; remove missing values
	if (any(ismissing(oneD_data_in)))
		ii = ind(.not. ismissing(oneD_data_in))  ; find locations of non-missing values
		oneD_data = oneD_data_in(ii)            ; new array containing no missing values
	else
		oneD_data = oneD_data_in 
	end if

	qsort(oneD_data)				; sort data - ascending order
	numel = dimsizes(oneD_data)		; find number of elements
	
	; declare stats array, of type data, to hold values
	stats = new(5, typeof(data))
	
	perc_point = (/ 5., 25., 50., 75., 95. /)
	
	do j = 0, 4
	
		; find point of percentile - perc_point must be between 0 and 100
		pt_x = round(.01*perc_point(j)*(numel-1),3)
		
		; quick correction, in case of negative pt_x
		if (pt_x .lt. 0)
			print("Catching error in stats calc due to small sample size")
			pt_x = 0
		end if
		
		stats(j) = oneD_data(pt_x)

	end do
	
	return(stats)
end

;-------------------------------------------------------------------
; Calc all stats:
; calculates, standard stats for domain averages5th, 25th , median, 75th and 95th percentiles. 
; Returns them an array of size 9:
;	stats(0) = mean(avg)
;	stats(1) = stdev
; 	stats(2) = min
;	stats(3) = max 
;	stats(4) = median
; 	stats(5) = lower quartile
;	stats(6) = upper quartile
;	stats(7) = 5th percentile
;	stats(8) = 95th percentile

; Input is a single array of data. The data may be multi- or single-dimensional. 

undef("calc_standard_stats")
function calc_standard_stats(data:numeric)
local oneD_data, ndims, numel, oneD_data_in, ii, boxplot_stats
begin

	ndims = dimsizes(dimsizes(data))
	if (ndims .gt. 1)
		oneD_data_in = ndtooned(data)
	else
		oneD_data_in = data 
	end if

	; declare stats array, of type data, to hold values
	stats = new(9, typeof(data), -9999)

	; remove missing values
	if (any(ismissing(oneD_data_in)))
		if (any(.not.ismissing(oneD_data_in)))
			ii = ind(.not. ismissing(oneD_data_in))  ; find locations of non-missing values
			oneD_data = oneD_data_in(ii)            ; new array containing no missing values
		else	; if no data then let's get out of here
			stats = 0.0
			return(stats)
		end if
	else
		oneD_data = oneD_data_in
	end if

	qsort(oneD_data)				; sort data - ascending order
	numel = dimsizes(oneD_data)		; find number of elements
	
	
	; calculate stats
	stats(0) = avg(oneD_data)    ; Average/mean
	stats(1) = stddev(oneD_data) ; standard deviation
	stats(2) = min(oneD_data)    ; minimum
	stats(3) = max(oneD_data)    ; maximum
	
	; Calculate boxplot stats, percentiles =  (5th, 25th, 50th, 75th, 95th) 
	boxplot_stats = calc_boxplot_stats(oneD_data)
	
	stats(4) = boxplot_stats(2)  ; median
	stats(5) = boxplot_stats(1)  ; lower quartile
	stats(6) = boxplot_stats(3)  ; upper quartile
	stats(7) = boxplot_stats(0)  ; 5th percentile
	stats(8) = boxplot_stats(4)  ; 95th percentile
	
	return(stats)
end

;-------------------------------------------------------------------

; !!!WARNING: the geometric stat calculations run the risk of dividing by 0,
; Caution should be used run running with.
; calculate geometric mean
; mu_g = (A_1*A_2*...*A_N)^(1/N) = exp( SUM(ln(A_i))/N )

undef("calc_geom_mean")
function calc_geom_mean(data:numeric)
local oneD_data, numel, log_dat_sq
begin

	sum_log = sum(log(data))		; sum log of data
	oneD_data = ndtooned(data)		; put data into 1 dimension
	numel = dimsizes(oneD_data)	; find number of elements

	mu_g = exp(sum_log/numel)
	return(mu_g)
	
end

; Calculate geometric standard deviation from array of data (sigma_g)
; sigma_g = exp(sqrt(SIGMA_i{ln(A_i/mu_g)}^2/N))

undef("calc_geom_stddev")
function calc_geom_stddev(data:numeric)
local mu_g, oneD_data, numel, log_dat_sq
begin
	mu_g = calc_geom_mean(data)			; calculate average
	oneD_data = ndtooned(data)
	numel = dimsizes(oneD_data)	; Find N
	
	if (mu_g .eq. 0.0)
		sigma_g = 0
	else
	
		log_dat_sq = log(data/mu_g) ^ 2 ; calculate bit that goes inside sum	
		sigma_g = exp( sqrt( sum(log_dat_sq)/numel) )
	end if
	
	return(sigma_g)
end

; Calculate geometric standard deviation from array of data (sigma_g)
; weighted by air mass. Needs air mass (AM) in same dimensions as data
; sigma_g = exp(sqrt(SIGMA_i{ln[weight*(A_i/mu_g)^2]/total_weighting))
; SAN - is this really meaningful?... - need to revise...

;undef("calc_geom_stddev_wt")
;function calc_geom_stddev_wt(data:numeric,AM:numeric)
;local mu_g, log_dat_sq,  tot_AM
;begin
;	tot_AM = sum(AM)			; find total air mass
;	mu_g = sum(data*AM)/tot_AM			; calculate weighted average
;	
;	log_dat_sq = log(data/mu_g) ^ 2 	; calculate bit that goes inside sum
;	
;	sigma_g_wt = exp( sqrt( sum(log_dat_sq*AM)/tot_AM) )
;	return(sigma_g_wt)
;end

;--------------------------------------------------------------------------------


; Obtain the unique values in an array (this is a stop-gap until we get upgrade
; to NCL 6.3.0 or higher)
;
undef("get_local_unique")
function get_local_unique(data:numeric)
local oneD_temp, oneD_data, data_store, numel, i, icount
begin
	; collapse the data array to 1D
	oneD_temp = ndtooned(data)
	; get rid of missing values
	oneD_data = rm_missing_1d(oneD_temp)

	numel = dimsizes(oneD_data)

	data_store = new(numel,typeof(data),data@_FillValue)

	icount = 0
	do i = 0, numel-1
		if(ismissing(oneD_data(i))) 
			continue
		end if
		if(.not.any(oneD_data(i) .eq. data_store))
			data_store(icount) = oneD_data(i)
			icount = icount + 1
		end if
	end do

	return(data_store(0:(icount-1)))

end




end