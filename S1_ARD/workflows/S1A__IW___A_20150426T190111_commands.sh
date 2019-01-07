# this script was created automatically by pyroSAR on Mon Jan 07 23:09:36 2019

export base=/geonfs01_vol1/ve39vem/test/dem_test/process/S1A_IW_GRDH_1SDV_20150426T190111_20150426T190136_005659_00741C_FFB7.SAFE

/cluster/GAMMA_SOFTWARE-20180704/ISP/bin/par_S1_GRD $base/measurement/s1a-iw-grd-vh-20150426t190111-20150426t190136-005659-00741c-002.tiff $base/annotation/s1a-iw-grd-vh-20150426t190111-20150426t190136-005659-00741c-002.xml $base/annotation/calibration/calibration-s1a-iw-grd-vh-20150426t190111-20150426t190136-005659-00741c-002.xml - $base/S1A__IW___A_20150426T190111_VH_grd.par $base/S1A__IW___A_20150426T190111_VH_grd - - - - -

/cluster/GAMMA_SOFTWARE-20180704/ISP/bin/par_S1_GRD $base/measurement/s1a-iw-grd-vv-20150426t190111-20150426t190136-005659-00741c-001.tiff $base/annotation/s1a-iw-grd-vv-20150426t190111-20150426t190136-005659-00741c-001.xml $base/annotation/calibration/calibration-s1a-iw-grd-vv-20150426t190111-20150426t190136-005659-00741c-001.xml - $base/S1A__IW___A_20150426T190111_VV_grd.par $base/S1A__IW___A_20150426T190111_VV_grd - - - - -

/cluster/GAMMA_SOFTWARE-20180704/ISP/bin/S1_OPOD_vec $base/S1A__IW___A_20150426T190111_VH_grd.par $base/osv/POEORB/S1A_OPER_AUX_POEORB_OPOD_20150517T123037_V20150425T225944_20150427T005944.EOF -

/cluster/GAMMA_SOFTWARE-20180704/ISP/bin/S1_OPOD_vec $base/S1A__IW___A_20150426T190111_VV_grd.par $base/osv/POEORB/S1A_OPER_AUX_POEORB_OPOD_20150517T123037_V20150425T225944_20150427T005944.EOF -

/cluster/GAMMA_SOFTWARE-20180704/ISP/bin/multi_look_MLI $base/S1A__IW___A_20150426T190111_VH_grd $base/S1A__IW___A_20150426T190111_VH_grd.par $base/S1A__IW___A_20150426T190111_VH_grd_mli $base/S1A__IW___A_20150426T190111_VH_grd_mli.par 11 9 - - -

/cluster/GAMMA_SOFTWARE-20180704/ISP/bin/multi_look_MLI $base/S1A__IW___A_20150426T190111_VV_grd $base/S1A__IW___A_20150426T190111_VV_grd.par $base/S1A__IW___A_20150426T190111_VV_grd_mli $base/S1A__IW___A_20150426T190111_VV_grd_mli.par 11 9 - - -

/cluster/GAMMA_SOFTWARE-20180704/DIFF/bin/gc_map $base/S1A__IW___A_20150426T190111_VH_grd_mli.par - /geonfs01_vol1/ve39vem/test/dem_test/demfile_gamma_SRTM-3Sec.par /geonfs01_vol1/ve39vem/test/dem_test/demfile_gamma_SRTM-3Sec $base/S1A__IW___A_20150426T190111_dem_seg.par $base/S1A__IW___A_20150426T190111_dem_seg $base/S1A__IW___A_20150426T190111_lut_coarse 1.0297146878 0.995627832068 - - - $base/S1A__IW___A_20150426T190111_inc - $base/S1A__IW___A_20150426T190111_pix $base/S1A__IW___A_20150426T190111_ls_map 8 0 -

/cluster/GAMMA_SOFTWARE-20180704/DIFF/bin/pixel_area $base/S1A__IW___A_20150426T190111_VH_grd_mli.par $base/S1A__IW___A_20150426T190111_dem_seg.par $base/S1A__IW___A_20150426T190111_dem_seg $base/S1A__IW___A_20150426T190111_lut_coarse $base/S1A__IW___A_20150426T190111_ls_map $base/S1A__IW___A_20150426T190111_inc $base/S1A__IW___A_20150426T190111_pixel_area_fine - - -

/cluster/GAMMA_SOFTWARE-20180704/ISP/bin/radcal_MLI $base/S1A__IW___A_20150426T190111_VH_grd_mli $base/S1A__IW___A_20150426T190111_VH_grd_mli.par - $base/S1A__IW___A_20150426T190111_VH_grd_mli_cal - - - 1 - - $base/S1A__IW___A_20150426T190111_ellipse_pixel_area

/cluster/GAMMA_SOFTWARE-20180704/LAT/bin/ratio $base/S1A__IW___A_20150426T190111_ellipse_pixel_area $base/S1A__IW___A_20150426T190111_pixel_area_fine $base/S1A__IW___A_20150426T190111_ratio_sigma0 2858 1 1 -

/cluster/GAMMA_SOFTWARE-20180704/LAT/bin/product $base/S1A__IW___A_20150426T190111_VH_grd_mli $base/S1A__IW___A_20150426T190111_ratio_sigma0 $base/S1A__IW___A_20150426T190111_VH_grd_mli_pan 2858 1 1 -

/cluster/GAMMA_SOFTWARE-20180704/DIFF/bin/geocode_back $base/S1A__IW___A_20150426T190111_VH_grd_mli_pan 2858 $base/S1A__IW___A_20150426T190111_lut_coarse $base/S1A__IW___A_20150426T190111_VH_grd_mli_pan_geo 3118 - 2 - - - -

/cluster/GAMMA_SOFTWARE-20180704/LAT/bin/lin_comb 1 $base/S1A__IW___A_20150426T190111_VH_grd_mli_pan_geo 0 0.776698728308 $base/S1A__IW___A_20150426T190111_VH_grd_mli_pan_geo_flat 3118 - - - - -

/cluster/GAMMA_SOFTWARE-20180704/LAT/bin/sigma2gamma $base/S1A__IW___A_20150426T190111_VH_grd_mli_pan_geo_flat $base/S1A__IW___A_20150426T190111_inc $base/S1A__IW___A_20150426T190111_VH_grd_mli_norm_geo 3118

/cluster/GAMMA_SOFTWARE-20180704/LAT/bin/product $base/S1A__IW___A_20150426T190111_VV_grd_mli $base/S1A__IW___A_20150426T190111_ratio_sigma0 $base/S1A__IW___A_20150426T190111_VV_grd_mli_pan 2858 1 1 -

/cluster/GAMMA_SOFTWARE-20180704/DIFF/bin/geocode_back $base/S1A__IW___A_20150426T190111_VV_grd_mli_pan 2858 $base/S1A__IW___A_20150426T190111_lut_coarse $base/S1A__IW___A_20150426T190111_VV_grd_mli_pan_geo 3118 - 2 - - - -

/cluster/GAMMA_SOFTWARE-20180704/LAT/bin/lin_comb 1 $base/S1A__IW___A_20150426T190111_VV_grd_mli_pan_geo 0 0.776698728308 $base/S1A__IW___A_20150426T190111_VV_grd_mli_pan_geo_flat 3118 - - - - -

/cluster/GAMMA_SOFTWARE-20180704/LAT/bin/sigma2gamma $base/S1A__IW___A_20150426T190111_VV_grd_mli_pan_geo_flat $base/S1A__IW___A_20150426T190111_inc $base/S1A__IW___A_20150426T190111_VV_grd_mli_norm_geo 3118

/cluster/GAMMA_SOFTWARE-20180704/LAT/bin/linear_to_dB $base/S1A__IW___A_20150426T190111_VH_grd_mli_norm_geo $base/S1A__IW___A_20150426T190111_VH_grd_mli_norm_geo_db 3118 0 -99

/cluster/GAMMA_SOFTWARE-20180704/DISP/bin/data2geotiff $base/S1A__IW___A_20150426T190111_dem_seg.par $base/S1A__IW___A_20150426T190111_VH_grd_mli_norm_geo_db 2 /geonfs01_vol1/ve39vem/test/dem_test/S1A__IW___A_20150426T190111_VH_grd_mli_norm_geo_db.tif -99

/cluster/GAMMA_SOFTWARE-20180704/LAT/bin/linear_to_dB $base/S1A__IW___A_20150426T190111_VV_grd_mli_norm_geo $base/S1A__IW___A_20150426T190111_VV_grd_mli_norm_geo_db 3118 0 -99

/cluster/GAMMA_SOFTWARE-20180704/DISP/bin/data2geotiff $base/S1A__IW___A_20150426T190111_dem_seg.par $base/S1A__IW___A_20150426T190111_VV_grd_mli_norm_geo_db 2 /geonfs01_vol1/ve39vem/test/dem_test/S1A__IW___A_20150426T190111_VV_grd_mli_norm_geo_db.tif -99

