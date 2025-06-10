from GBT_waterfall import * 
from dysh.fits.gbtfitsload import GBTFITSLoad, GBTOffline
import warnings
from astropy.utils.exceptions import AstropyWarning

import logging, sys
logging.disable(sys.maxsize)

if __name__ == "__main__":
    

    PF800_RFI = ["TRFI_010125_81",
                        "TRFI_011625_81",
                        "TRFI_021425_81",
                        "TRFI_022125_81",
                        "TRFI_050424_81",
                        "TRFI_050524_81",
                        "TRFI_050624_81",
                        "TRFI_050924_81",
                        "TRFI_051124_81",
                        "TRFI_051224_81",
                        "TRFI_052724_81",
                        "TRFI_053124_81",
                        "TRFI_101924_81",
                        "TRFI_112924_81",
                        "TRFI_121324_81",
                        "TRFI_122124_81",
                        "TRFI_122624_81",
                        "TRFI_122824_81",
                        "TRFI_123124_81"]

    outdir_raw = "/users/jfrothin/RFI/plots/raw/"
    outdir_med = "/users/jfrothin/RFI/plots/median/"

    # for session_ID in test_list:
    #     sdf = GBTOffline(session_ID) # for data in /home/sdfits
        
    #     GBT_waterfall(sdf, session_ID, 
    #                   fmin_GHz=1.97, 
    #                   fmax_GHz=2, 
    #                   band_allocation="starlink", 
    #                   outdir=outdir, 
    #                   cal_type="raw_data")
        
    #     GBT_waterfall(sdf, session_ID, 
    #                   fmin_GHz=1.97, 
    #                   fmax_GHz=2, 
    #                   band_allocation="starlink", 
    #                   outdir=outdir, 
    #                   cal_type="median_subtract")
    
    session_ID = PF800_RFI[0]
    sdf = GBTOffline(session_ID) # for data in /home/sdfits
    
    GBT_waterfall(sdf, session_ID,  
                  fmin_GHz=0.75,
                  outdir=outdir_raw, 
                  band_allocation="25B-184",
                  cal_type="raw_data")
    
    GBT_waterfall(sdf, session_ID,
                  fmin_GHz=0.75,
                  outdir=outdir_med, 
                  band_allocation="25B-184",
                  cal_type="median_subtract")