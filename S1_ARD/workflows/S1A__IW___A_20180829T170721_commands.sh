# this script was created automatically by pyroSAR on Tue Mar 12 14:27:20 2019

export base=/home/truc_jh/Desktop/S1_ARD/data/GAMMA/stack/process/S1A_IW_GRDH_1SDV_20180829T170721_20180829T170746_023464_028DE0_5310.SAFE

par_S1_GRD $base/measurement/s1a-iw-grd-vh-20180829t170721-20180829t170746-023464-028de0-002.tiff $base/annotation/s1a-iw-grd-vh-20180829t170721-20180829t170746-023464-028de0-002.xml $base/annotation/calibration/calibration-s1a-iw-grd-vh-20180829t170721-20180829t170746-023464-028de0-002.xml - $base/S1A__IW___A_20180829T170721_VH_grd.par $base/S1A__IW___A_20180829T170721_VH_grd - - - - -

par_S1_GRD $base/measurement/s1a-iw-grd-vv-20180829t170721-20180829t170746-023464-028de0-001.tiff $base/annotation/s1a-iw-grd-vv-20180829t170721-20180829t170746-023464-028de0-001.xml $base/annotation/calibration/calibration-s1a-iw-grd-vv-20180829t170721-20180829t170746-023464-028de0-001.xml - $base/S1A__IW___A_20180829T170721_VV_grd.par $base/S1A__IW___A_20180829T170721_VV_grd - - - - -

S1_OPOD_vec $base/S1A__IW___A_20180829T170721_VH_grd.par $base/osv/POEORB/S1A_OPER_AUX_POEORB_OPOD_20180918T120731_V20180828T225942_20180830T005942.EOF -

S1_OPOD_vec $base/S1A__IW___A_20180829T170721_VV_grd.par $base/osv/POEORB/S1A_OPER_AUX_POEORB_OPOD_20180918T120731_V20180828T225942_20180830T005942.EOF -

multi_look_MLI $base/S1A__IW___A_20180829T170721_VH_grd $base/S1A__IW___A_20180829T170721_VH_grd.par $base/S1A__IW___A_20180829T170721_VH_grd_mli $base/S1A__IW___A_20180829T170721_VH_grd_mli.par 2 2 - - -

multi_look_MLI $base/S1A__IW___A_20180829T170721_VV_grd $base/S1A__IW___A_20180829T170721_VV_grd.par $base/S1A__IW___A_20180829T170721_VV_grd_mli $base/S1A__IW___A_20180829T170721_VV_grd_mli.par 2 2 - - -

gc_map $base/S1A__IW___A_20180829T170721_VH_grd_mli.par - /home/truc_jh/Desktop/S1_ARD/data/DEM/alps_dem_gamma_AW3D30.par /home/truc_jh/Desktop/S1_ARD/data/DEM/alps_dem_gamma_AW3D30 $base/S1A__IW___A_20180829T170721_dem_seg_geo.par $base/S1A__IW___A_20180829T170721_dem_seg_geo $base/S1A__IW___A_20180829T170721_lut_init 1.0 1.0 - - - $base/S1A__IW___A_20180829T170721_inc_geo - $base/S1A__IW___A_20180829T170721_pix_geo $base/S1A__IW___A_20180829T170721_ls_map_geo 8 2 -

pixel_area $base/S1A__IW___A_20180829T170721_VH_grd_mli.par $base/S1A__IW___A_20180829T170721_dem_seg_geo.par $base/S1A__IW___A_20180829T170721_dem_seg_geo $base/S1A__IW___A_20180829T170721_lut_init $base/S1A__IW___A_20180829T170721_ls_map_geo $base/S1A__IW___A_20180829T170721_inc_geo - - - - $base/S1A__IW___A_20180829T170721_pix_fine -

product $base/S1A__IW___A_20180829T170721_VH_grd_mli $base/S1A__IW___A_20180829T170721_pix_fine $base/S1A__IW___A_20180829T170721_VH_grd_mli_pan 16101 1 1 -

geocode_back $base/S1A__IW___A_20180829T170721_VH_grd_mli_pan 16101 $base/S1A__IW___A_20180829T170721_lut_init $base/S1A__IW___A_20180829T170721_VH_grd_mli_pan_geo 5927 - 2 - - - -

sigma2gamma $base/S1A__IW___A_20180829T170721_VH_grd_mli_pan_geo $base/S1A__IW___A_20180829T170721_inc_geo $base/S1A__IW___A_20180829T170721_VH_grd_mli_norm_geo 5927

product $base/S1A__IW___A_20180829T170721_VV_grd_mli $base/S1A__IW___A_20180829T170721_pix_fine $base/S1A__IW___A_20180829T170721_VV_grd_mli_pan 16101 1 1 -

geocode_back $base/S1A__IW___A_20180829T170721_VV_grd_mli_pan 16101 $base/S1A__IW___A_20180829T170721_lut_init $base/S1A__IW___A_20180829T170721_VV_grd_mli_pan_geo 5927 - 2 - - - -

sigma2gamma $base/S1A__IW___A_20180829T170721_VV_grd_mli_pan_geo $base/S1A__IW___A_20180829T170721_inc_geo $base/S1A__IW___A_20180829T170721_VV_grd_mli_norm_geo 5927

linear_to_dB $base/S1A__IW___A_20180829T170721_VH_grd_mli_norm_geo $base/S1A__IW___A_20180829T170721_VH_grd_mli_norm_geo_db 5927 0 -99

data2geotiff $base/S1A__IW___A_20180829T170721_dem_seg_geo.par $base/S1A__IW___A_20180829T170721_VH_grd_mli_norm_geo_db 2 /home/truc_jh/Desktop/S1_ARD/data/GAMMA/stack/S1A__IW___A_20180829T170721_VH_grd_mli_norm_geo_db.tif -99

linear_to_dB $base/S1A__IW___A_20180829T170721_VV_grd_mli_norm_geo $base/S1A__IW___A_20180829T170721_VV_grd_mli_norm_geo_db 5927 0 -99

data2geotiff $base/S1A__IW___A_20180829T170721_dem_seg_geo.par $base/S1A__IW___A_20180829T170721_VV_grd_mli_norm_geo_db 2 /home/truc_jh/Desktop/S1_ARD/data/GAMMA/stack/S1A__IW___A_20180829T170721_VV_grd_mli_norm_geo_db.tif -99

data2geotiff $base/S1A__IW___A_20180829T170721_dem_seg_geo.par $base/S1A__IW___A_20180829T170721_inc_geo 2 /home/truc_jh/Desktop/S1_ARD/data/GAMMA/stack/S1A__IW___A_20180829T170721_inc_geo.tif -99

data2geotiff $base/S1A__IW___A_20180829T170721_dem_seg_geo.par $base/S1A__IW___A_20180829T170721_ls_map_geo 5 /home/truc_jh/Desktop/S1_ARD/data/GAMMA/stack/S1A__IW___A_20180829T170721_ls_map_geo.tif 0
