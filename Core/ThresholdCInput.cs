using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using uQlustCore.Interface;
using uQlustCore.Profiles;

namespace uQlustCore
{
    public class ThresholdCInput:BaseCInput,IAutomaticProfiles
    {
        [Description("Threshold clustering: What is the minimal number of structures in the cluster")]
        public int bakerNumberofStruct;
        [Description("Threshold clustering algorithm: set the value of distance threshold, below this value structure are treated as similar")]
        public float distThresh;
        [Description("Distance measures used for threshold clustering")]
        public DistanceMeasures hDistance = DistanceMeasures.HAMMING;
        [Description("Profiles used for threshold custering when weihgted hamming is used")]
        public string hammingProfile;
        [Description("In case of RMSD or MAXSUB which atoms should be used in threshold clustering")]
        public PDB.PDBMODE hAtoms = PDB.PDBMODE.ONLY_CA;
        [Description("Use 1DJury to find reference vectors")]
        public bool reference1Djury;

        public void GenerateAutomaticProfiles(string fileName)
        {
            ProfileTree t = ProfileAutomatic.AnalyseProfileFile(fileName, SIMDIST.DISTANCE);
            string profileName = "automatic_distance.profile";
            t.SaveProfiles(profileName);
            hammingProfile = profileName;           
        }

    }
}
