#include "TH1.h"
#include "TString.h"

#include "RooWorkspace.h"

RooWorkspace* getWorkspace(TH1* ww,
                           TH1* top,
                           TH1* dytt,
                           TH1* vv, 
                           TH1* sig,
                           TString label);
