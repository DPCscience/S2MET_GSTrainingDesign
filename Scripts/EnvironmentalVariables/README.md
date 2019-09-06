## Notes on Collecting Environmental Data

This database returns variables per horizon per component (these can be found in a description [here](https://sdmdataaccess.nrcs.usda.gov/documents/TableColumnDescriptionsReport.pdf)). The following variables will be kept:


  Variable      | Definition
----------------|---------------------------------------------------------------
hzdept_r        | Distance from top of soil to top of horizon (cm)
hzdepb_r        | Distance from top of soil to bottom of horizon (cm)
sandtotal_r     | Mineral particles 0.05mm to 2.0mm in equivalent diameter as a weight 
                | percentage of the less than 2 mm fraction.
silttotal_r     | Mineral particles 0.002 to 0.05mm in equivalent diameter as a weight 
                | percentage of the less than 2.0mm fraction.
claytotal_r     | Mineral particles less than 0.002mm in equivalent diameter as a weight 
                | percentage of the less than 2.0mm fraction.
om_r            | The amount by weight of organic matter expressed as a weight percentage 
                | of the less than 2 mm soil material
ph1to1h2o_r     | pH


Query the National Soil Database for site information in Canada

The analagous variables from the USDA Soil Survey in the National Soil DataBase (NSDB, Canada) are:

  Variable      | NSDB Analogue
----------------|---------------------------------------------------------------
hzdept_r        | UDEPTH
hzdepb_r        | LDEPTH - UDEPTH
sandtotal_r     | TSAND
silttotal_r     | TSILT
claytotal_r     | TCLAY
om_r            | ORGCARB
ph1to1h2o_r     | PH2

