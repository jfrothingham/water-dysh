from GBT_waterfall import * 
from dysh.fits.gbtfitsload import GBTFITSLoad, GBTOffline

if __name__ == "__main__":
    # uwbr
    ryan_uwbr_tests = ["TGBT25A_607_01", "TGBT25A_607_02", "TGBT25A_607_03", "TGBT25A_607_04",]
    uwbr_rfi_scans = ["TRFI_050825_251", "TRFI_051524_251", "TRFI_111424_251"]
    ods_tests = [#"TGBT25A_607_04", # didn't plot for some reason
                 #"TGBT25A_607_05", # didn't plot for some reason
                 "TGBT25A_607_06", 
                 "TGBT25A_607_07", 
                 "TGBT25A_607_08",
                 "TGBT25A_607_09",
                 "TGBT25A_607_10",
                 "TGBT25A_607_11",
                #  "TGBT25A_607_12", # ran in vegas pulsar mode
                #  "TGBT25A_606_08", # cyclic spectroscopy (pulsar mode data)
                #  "TGBT25A_606_09", # cyclic spectroscopy (pulsar mode data)
                #  "TGBT25A_606_10", # cyclic spectroscopy (pulsar mode data)
                #  "TGBT25A_606_11", # cyclic spectroscopy (pulsar mode data)
                #  "TGBT25A_606_12", # cyclic spectroscopy (pulsar mode data)
                 ]

    # x-band
    remijan_rfi = ["AGBT25A_218_23"]
    nrdz_spaceX_tests = ["TGBT24B_609_01","TGBT24B_609_02","TGBT24B_609_03","TGBT24B_609_05","TGBT24B_609_05"]
    xband_rfi_scans = ["TRFI_020325_X1", # completed -- do not need to run again
            #"TRFI_033025_X1", # this one broke -- look into why
            # "TRFI_040625_X1", # this one broke -- look into why
            "TRFI_041125_X1",
            "TRFI_050424_X1",
            "TRFI_050524_X1",
            "TRFI_050624_X1",
            "TRFI_050924_X1",
            "TRFI_051124_X1",
            "TRFI_051224_X1",
            "TRFI_051824_X1",
            "TRFI_060124_X1",
            "TRFI_060224_X1",
            "TRFI_060224_X2",
            "TRFI_060324_X1",
            "TRFI_101924_X1",
            "TRFI_112524_X1",
            "TRFI_120624_X1",
            "TRFI_122824_X1",
            # "TRFI_122924_X1", # this one broke -- look into why
            "TRFI_123124_X1"
            ]
    older_xband = ["TRFI_123124_X1"]

    # ku-band 
    ku_band_post_spaceX_sessions = ["TRFI_123124_U1", "TRFI_122424_U1", "TRFI_122224_U1","TRFI_122024_U1"]
    
    from_archive = False
    outdir = "/home/scratch/dbautist/operational_data_sharing/plots/waterfalls/"

    for session_ID in ods_tests:
        sdf = GBTOffline(session_ID) # for data in /home/sdfits
        
        GBT_waterfall(sdf, session_ID, 
                      fmin_GHz=1.97, 
                      fmax_GHz=2, 
                      band_allocation="starlink", 
                      outdir=outdir, 
                      cal_type="raw_data")
        
        GBT_waterfall(sdf, session_ID, 
                      fmin_GHz=1.97, 
                      fmax_GHz=2, 
                      band_allocation="starlink", 
                      outdir=outdir, 
                      cal_type="median_subtract")
