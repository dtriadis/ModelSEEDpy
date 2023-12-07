# -*- coding: utf-8 -*-

from __future__ import absolute_import

# set the warning format to be on a single line
import sys
import logging
import cobra
import warnings as _warnings
from os import name as _name
from os.path import abspath as _abspath
from os.path import dirname as _dirname
from modelseedpy_freiburgermsu.helpers import config

__author__ = "Christopher Henry"
__email__ = "chenry@anl.gov"
__version__ = "0.2.2"

logger = logging.getLogger(__name__)

print("modelseedpy_freiburgermsu", __version__)

if sys.version_info[0] == 2:
    logger.warning(
        "Python 2 is reaching end of life (see "
        "https://www.python.org/dev/peps/pep-0373/) and many cobra "
        "dependencies have already dropped support. At the moment it *should* "
        "still work but we will no longer actively maintain Python 2 support."
    )

if "e0" not in cobra.medium.annotations.compartment_shortlist["e"]:
    cobra.medium.annotations.compartment_shortlist["e"].append("e0")

import modelseedpy_freiburgermsu
from modelseedpy_freiburgermsu.core import (
    RastClient, MSGenome, MSBuilder, MSMedia, MSGrowthPhenotypes, MSModelUtil,
    FBAHelper, MSEditorAPI, MSATPCorrection, MSGapfill, MSEquation, OptlangHelper,
    commscores_report, steadycom_report
)
from modelseedpy.core.exceptions import *

from modelseedpy.community import (
    MSCommunity, CommunityMember, MSKineticsFBA, CommScores, CommPhitting,
    MSSteadyCom, build_from_species_models, phenotypes)

from modelseedpy_freiburgermsu.biochem import ModelSEEDBiochem

from modelseedpy_freiburgermsu.fbapkg import (
    BaseFBAPkg, RevBinPkg, ReactionUsePkg, SimpleThermoPkg, TotalFluxPkg, BilevelPkg, CommKineticPkg,
    KBaseMediaPkg, FluxFittingPkg, ProteomeFittingPkg, GapfillingPkg, MetaboFBAPkg, FlexibleBiomassPkg,
    ProblemReplicationPkg, FullThermoPkg, MSPackageManager, ObjConstPkg, ChangeOptPkg, ElementUptakePkg
)

from modelseedpy_freiburgermsu.multiomics import MSExpression
