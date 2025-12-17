"""This script create a study, which mean a set of repository
which can be convenient to manage path. (might doesn't work on windows...)
"""
from brainpipe.system import study

if __name__ == "__main__":
    # Add parent directory to path to import config
    import sys
    from pathlib import Path
    current_dir = Path(__file__).resolve().parent
    # olfaction_scripts_old is 2 levels up from 0_scripts_Etienne
    root_dir = current_dir.parent.parent
    sys.path.append(str(root_dir))
    
    try:
        import config
        path = config.get_study_path()
        name = config.STUDY_NAME
        print(f"-> Using Data Path from config: {path}")
    except ImportError:
        print("-> Config not found, falling back to local relative path")
        # Fallback if config isn't found
        path = str(current_dir.parent) # 3_BrainpipeScripts
        name = 'Olfacto'

    # Create the study :
    st = study(name)
    st.add(path)
