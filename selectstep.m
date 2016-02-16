function stepoptions = selectstep(step)

% selectstep returns a structure containing the appropriate step options
% related to the step number chosen

% Init
CHOOSE_RUN = [];
CHECK_SYNC = [];
SHOW_WFbefAnalysis = [];
SHOW_allWF = [];
SHOW_1WFperfile = [];
SHOW_WFref = [];

switch step
    case '1'
        CHOOSE_RUN = 1;
    case '2'
        CHOOSE_RUN = 0;
        CHECK_SYNC = 1;
    case '3'
        CHOOSE_RUN = 0;
        CHECK_SYNC = 0;
        SHOW_WFbefAnalysis = 1;
    case '4'
        CHOOSE_RUN = 0;
        CHECK_SYNC = 0;
        SHOW_WFbefAnalysis = 0;
        SHOW_WFref = 1;
    case '5'
        CHOOSE_RUN = 0;
        CHECK_SYNC = 0;
        SHOW_WFbefAnalysis = 0;
        SHOW_allWF = 1;
        SHOW_1WFperfile = 0;
        SHOW_WFref = 0;
    case '6'
        CHOOSE_RUN = 0;
        CHECK_SYNC = 0;
        SHOW_WFbefAnalysis = 0;
        SHOW_allWF = 0;
        SHOW_1WFperfile = 1;
        SHOW_WFref = 0;
end

% create structure
stepoptions.CHOOSE_RUN = CHOOSE_RUN;
stepoptions.CHECK_SYNC = CHECK_SYNC;
stepoptions.SHOW_allWF = SHOW_allWF;
stepoptions.SHOW_1WFperfile = SHOW_1WFperfile;
stepoptions.SHOW_WFref = SHOW_WFref;
stepoptions.SHOW_WFbefAnalysis = SHOW_WFbefAnalysis;

end



