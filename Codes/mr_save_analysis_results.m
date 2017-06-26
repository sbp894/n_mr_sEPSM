function mr_save_analysis_results(PowerModCell,resultsDir,resultPostfix,paramsIN)


Library.parsave([resultsDir 'psd' filesep 'psd' resultPostfix '.mat'], PowerModCell);
    
Library.parsave([resultsDir 'paramsIN' filesep 'paramsIN' resultPostfix '.mat'], paramsIN);