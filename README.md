# Olfaction-sEEG-behav projec
This is breakdown of old files: Data and Scripts.

olfaction_scripts_old/
│
├── 0_RawData_trc_avi_behav/                         # Raw Data
│   ├── S0X_iden_ordre=FM                            # Subject '0X' with identifier iden and data presented in order FM
│   │    ├── S0X_iden_Behavior                       # Subject '0X' Behav Data (configPos, header, ParaOlfacto, raw, res, xlsx files)
│   │    └── S0X_iden_SignalsSEEG                    # Subject '0X' sEEG Data (TRC files)
│   └── ...                                          # 6 Subjects in total
│
├── 1bis_OE_BaseSam/                                 # Behav results, figures and diverse scripts
│   ├── JPlailly201306_seeg_ALS                      # old scripts
│   │    ├── behavior
│   │    │    ├── plots                              # List of plots (png/pdf)
│   │    │    ├── respiration_amplitude              # numpy arrays (.npy files) 
│   │    │    ├── Compute_Plot_respi.ipynb           # plotting notebook
│   │    │    ├── Name.py                            # multiple python files for different anaylsis
│   │    │    └── ...                                # multiple python csv and xls files with results
│   │    ├── conn_analysis                   
│   │    │    └── Name.py                            # multiple python files for different anaylsis
│   │    ├── Scripts_concat
│   │    │    ├── Concatenate_sigs_Odor_Rest.ipynb   # notebook
│   │    │    ├── Export_correct_sig.ipynb           # notebook
│   │    │    └── ...                                # .m files
│   │    ├── seeg_analysis
│   │    │    ├── bytrial_envelop_cycles             # figures (png)
│   │    │    ├── diff_cyles_maps_0-40               # figures (png)
│   │    │    ├── Name.py                            # multiple python files for different anaylsis
│   │    │    └── ...                                # png, pdf and xls files
│   │    └── database JPlailly201306_seeg_ALS.sqlite # Not sure what is this
│   ├── OLD                         # old scripts (JPlailly201306_seeg_david and JPlailly201306_seeg_win)
│   └── Results_Behavior.py         # List of figures (Pdf, png files)
│        ├── Figures                # Subject '0X' Behav Data (configPos, header, ParaOlfacto, raw, res, xlsx files)
│        └── ...                    # xls and ods files: dataframes of results