import moments
import numpy as np
import os
#import Misc2 #this is an edited Misc frm moments with a change in bootstrap



# Input
VCFInput = "../yvb_gadma.recode.vcf"
pop_id = ["pop1"]
sampleSizes = [146]
save_dir = "boot"

# Make data dict
dd = moments.Misc.make_data_dict_vcf(VCFInput, "../pop1_map2.txt")

# Create output dir if it doesn't exist
os.makedirs(save_dir, exist_ok=True)

# Run bootstrap
moments.Misc.bootstrap(dd, pop_ids=pop_id, projections=sampleSizes,
                mask_corners=True, polarized=False, bed_filename=None,
                num_boots=101, save_dir=save_dir)

