using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using uQlustCore.Interface;
using uQlustCore.Profiles;

namespace uQlustCore
{
    public class RankingCInput:BaseCInput,IAutomaticProfiles
    {       
        public string alignFileName;
        public bool alignGenerate=false;
        [Description("Distance measure [OTHER]")]
        public DistanceMeasures oDistance = DistanceMeasures.HAMMING;
        public PDB.PDBMODE oAtoms = PDB.PDBMODE.ONLY_CA;
        [Description("Path to the profile for weighted hamming distance [OTHER]")]
        public string hammingProfile;
        [Description("Path to the profile for 1Djury [OTHER]")]
        public string juryProfile;
        [Description("Use 1DJury to find reference vectors")]
        public bool reference1Djury;
        [Description("Path to reference profile [OTHER]")]
        public string referenceProfile;

        public void GenerateAutomaticProfiles(string fileName)
        {
            ProfileTree t = ProfileAutomatic.AnalyseProfileFile(fileName, SIMDIST.DISTANCE);
            string profileName = "automatic_distance.profile";
            t.SaveProfiles(profileName);
            hammingProfile = profileName;
            t = ProfileAutomatic.AnalyseProfileFile(fileName, SIMDIST.SIMILARITY);
            profileName = "automatic_similarity.profile";
            t.SaveProfiles(profileName);
            juryProfile= profileName;
        }

    }
}
