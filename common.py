DIR_DATA = "/Users/lp2/data/projects/wildsupp"
DIR_AFS = "%s/seq/Joint-call/afs"%DIR_DATA
DIR_META = "%s/paper/meta"%DIR_DATA
DIR_TABLE = "%s/paper/tables"%DIR_DATA
DIR_PLOT = "%s/paper/plots"%DIR_DATA
DIR_QTL = "%s/paper/tables/qtls"%DIR_DATA
DIR_QTLPLOT = "%s/paper/plots/qtls"%DIR_DATA

SQTL_VERSION = "0.1"
STR_HEADER_CALL_REGIONS = """#sQTL version %s - call QTL regions output
#==========================================
# out_file=%s
# vcf_file=%s
# test sample=%s
# control sample=%s
# af_lenient=%.2f
# sd_lenient=%.1f
# af_stringent=%.2f
# sd_stringent=%.1f
# length_cutoff=%d
#==========================================
"""

import logging
import cPickle

logging.basicConfig(level=logging.DEBUG, format='%(asctime)s %(name)-12s: %(levelname)-8s %(message)s')
LOG = logging.getLogger("Wild suppressors")
