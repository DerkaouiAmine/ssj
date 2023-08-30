package fichier_Sob;
// methodes used :

/*
 * M0:CBC extend 3000
 * M1:CBC extend7000
 * M2:change only one vecteur per 7000 run and leave the others randoms
 * M3:construct the matrices element by element
 * M4:NAive Sob
 * M5:NAive Sob+Lms
 * M6:NAive sob+LMS 100 times
 * 
 * 
 * */
//PArametre for Best Wafom matrice=1;Worst Wafom matrice=0

public class WafomsStorage {
    private int[][] tabl;

    public WafomsStorage(int s, int methode,int BesOrWo) {
        if (s == 3) {
        	if ( methode == 0 && BesOrWo == 1) {
        		meilleursWafomsS3M0();        	}
        	if ( methode == 0 && BesOrWo == 0) {
        		worstWafomsS3M0();
        	}
        	
        	if ( methode == 1 && BesOrWo == 1) {
        		meilleursWafomsS3M1();        	}
        	if ( methode == 1 && BesOrWo == 0) {
        		worstWafomsS3M1();
        	}
        	
        	if ( methode == 2 && BesOrWo == 1) {
        		meilleursWafomsS3M2();        	}
        	if ( methode == 2 && BesOrWo == 0) {
        		worstWafomsS3M2();
        	}
        	
        	if ( methode == 3 && BesOrWo == 1) {
        		meilleursWafomsS3M3()  ;      	}
        	if ( methode == 3 && BesOrWo == 0) {
        		worstWafomsS3M3();
        	}
        	
        	if ( methode == 4 && BesOrWo == 1) {
        		meilleursWafomsS3M4() ;       	}
        	if ( methode == 4 && BesOrWo == 0) {
        		worstWafomsS3M4();
        	}
        	
        	
        	if ( methode == 5 && BesOrWo == 1) {
        		meilleursWafomsS3M5()   ;     	}
        	if ( methode == 5 && BesOrWo == 0) {
        		worstWafomsS3M5();
        	}
        	
        	if ( methode == 6 && BesOrWo == 1) {
        		meilleursWafomsS3M6() ;       	}
        	if ( methode == 6 && BesOrWo == 0) {
        		worstWafomsS3M6();
        	}
        	
        	
        	
        	
        	
            
       
        
        
        } else if (s == 6) {
           // meilleursWafomsS6();
        } // Ajoutez d'autres cas pour d'autres valeurs de s si nécessaire
    }
    
    
    /*
     * M0:CBC extend 3000
 
     * */
    private void meilleursWafomsS3M0() {
        
    	tabl = new int[][] {
    		{1046994105,532522000,152816856,75419729,53967615,19091783,15938120,4789585,2819841,1672586,623710,328735,139167,79238,39988,21016,14245,7347,2696,1508},
    		{1069297368,558378365,802900350,809188620,938537472,763178122,1042868645,940424680,662213238,676211638,832943665,942726162,1014808367,804978112,1069252715,660898319,914189279,807219605,584650230,754793274},
    		{1065712118,519591766,630359917,844464685,496438396,984524107,676074352,485055493,930891272,592170580,448320467,647096772,999816289,281575487,682941774,953289808,536779220,1043762743,624146450,456579996}

        };
    	
    }

    private void worstWafomsS3M0() {
    	tabl = new int[][] {
    		{544048359,284937534,260355067,98215387,40620890,20488884,12767131,8352775,2111812,1314218,1002770,389315,255195,124029,43128,18075,15021,4297,3852,1357},
    		{554621823,563548478,1028082886,770896851,537003194,806824459,756100253,913554581,759119390,668156570,776941497,560228320,942247245,837478900,738236909,907127963,776051840,910697921,1051098877,1054155163},
    		{537385419,281687899,1061718654,635020800,305091122,854930272,986635225,496411081,983115615,926570106,401426941,628011486,650653390,437513185,1045332176,785451650,384755544,584866170,918712457,472352341}
            
        };    }
    /*
   
     * M1:CBC extend7000

     * 
     * */

    private void meilleursWafomsS3M1() {
    	tabl = new int[][] {
    		{1046994105,398307896,181181815,95054231,34220517,21484503,12394050,4574640,2364507,1049891,926899,518558,221597,94087,34218,30883,14575,4451,2382,1950,554},
    		{1069297368,735226078,767106076,843169012,899961314,941641228,602368976,667437869,788484991,955352261,742762285,611944903,919851144,779148930,809466920,789232013,584593252,1051409679,983178797,677711189,996379055},
    		{1065712118,528293238,620731993,863648038,471982456,797987891,577453692,437494779,1031238497,885621485,347984506,622078176,850186234,393980218,972996136,1044507081,465061887,755970565,776187327,289374307,845703665}

            
        };    }

    private void worstWafomsS3M1() {
      
    	tabl = new int[][] {

{544048359,272122765,239806763,100730380,60707987,32541601,16216806,5332001,2255690,1673357,977125,380469,148902,89981,34166,28651,11995,5669,2470,1603,1007},
{554621823,545713435,1067514286,766073899,725571323,845208132,945632384,904144371,601641574,785359955,832528254,785251670,1044514597,738960223,1038704840,887364386,951668706,972058797,903831272,604955815,962139987},
{537385419,281480881,1063592985,611504236,526282491,857428469,799345022,470007140,574340618,1030188716,331112673,1031093835,627269216,347886913,764991417,829136862,498802332,671346433,911861371,345661266,568892714}
            
        };
    }
    /*
   
     * M2:change only one vecteur per 7000 run and leave the others randoms

     * 
     * */
    
    private void meilleursWafomsS3M2() {
    	tabl = new int[][] {
    		{1056608413,530967972,267738758,130270757,64876601,31948583,15157521,7108341,2211970,1071097,957737,470517,213152,124839,40353,26349,9908,6236,4051},
    		{912459566,971305106,904141872,556101454,874103430,1048728493,906778806,704838433,843936465,628551508,831951728,674551335,649175194,712302470,658038973,1035763900,944880364,944283062,941041629},
    		{657975599,344506383,945786601,721709468,314648839,641010118,850256436,411340742,782315336,890131885,293536464,915750056,800454010,535731598,1007137995,692201741,269843439,785937618,1071501392}

        };   }

    private void worstWafomsS3M2() {
    	tabl = new int[][] {
    		{1056608413,530967972,267738758,130270757,64876601,31948583,15157521,7108341,2211970,1071097,957737,470517,213152,124839,40353,26349,9908,6236,4051},
    		{912459566,971305106,904141872,556101454,874103430,1048728493,906778806,704838433,843936465,628551508,831951728,674551335,649175194,712302470,658038973,1035763900,944880364,944283062,941041629},
    		{657975599,344506383,945786601,721709468,314648839,641010118,850256436,411340742,782315336,890131885,293536464,915750056,800454010,535731598,1007137995,692201741,269843439,785937618,1071501392},

            
        };    }
    /*
 
     * M3:construct the matrices element by element

     * for this methodes we construct the matrice optimally so any other matrices is worst then the better one 
     * 
     * */
    private void meilleursWafomsS3M3() {
    	tabl = new int[][] {
    		{1070938843, 535753106, 267738343, 134188649, 65237965, 33522641, 15854459, 6790573, 4038514, 1418655, 974131, 277320, 199212, 117132, 40699, 27844, 11399, 4915, 3719, 2041},
    		{1058188086, 552349449, 1056987673, 768811916, 902998753, 992825677, 783694916, 735273998, 923051277, 816134121, 554761737, 772372081, 566657365, 911677216, 1032002328, 897874876, 1055045780, 923529241, 727348786, 622268637},
    		{1067074533, 511186581, 542757137, 1002221156, 347286984, 1011028352, 538502094, 497500292, 636366666, 927424421, 460901690, 1035843731, 632003044, 277321193, 956904058, 649641527, 307901518, 1065560734, 989526869, 430527640},
        };    }

    private void worstWafomsS3M3() {
    	System.out.println("for this methodes we construct the matrice optimally so any other matrices is worst then the better one");   }
    /*
     * M4:NAive Sob
    
     * 
     * */
    
    private void meilleursWafomsS3M4() {
    	tabl = new int[][] {
    		{536870912, 268435456, 134217728, 67108864, 33554432, 16777216, 8388608, 4194304, 2097152, 1048576, 524288, 262144, 131072, 65536, 32768, 16384, 8192, 4096, 2048, 1024, 512},
    		{536870912, 805306368, 671088640, 1006632960, 570425344, 855638016, 713031680, 1069547520, 538968064, 808452096, 673710080, 1010565120, 572653568, 858980352, 715816960, 1073725440, 536879104, 805318656, 671098880, 1006648320, 570434048},
    		{536870912, 268435456, 939524096, 738197504, 436207616, 1023410176, 562036736, 331350016, 975175680, 756023296, 431489024, 1072431104, 540672000, 271384576, 941195264, 742375424, 438312960, 1024462848, 565721088, 334244864, 976886272}
            
        };   }

    private void worstWafomsS3M4() {
    	System.out.println("for this methodes we construct the best one naively so we only take the result matrice as the best one");   }
	
    

    /*

     * M5:NAive Sob+Lms
     
     * 
     * */
    private void meilleursWafomsS3M5() {
    	tabl = new int[][] {
    	{615863594, 456192572, 222621575, 122457175, 46765358, 32698090, 11280134, 6295359, 2891556, 1113457, 795579, 297531, 174965, 108372, 48639, 22631, 13926, 7473, 2543, 1986, 690},
    	{1034796896, 593461357, 852653419, 676331799, 1057288269, 537357049, 817515700, 731203419, 1032013090, 593550387, 854993367, 677482897, 1059668788, 540701782, 815659376, 732331210, 1034787210, 593474869, 852646102, 676346834, 1057281278},
    	{758603405, 386556942, 829472416, 598622780, 530756666, 897059371, 739868085, 361307807, 850690919, 577284584, 492800743, 906923098, 755163853, 388137878, 828788190, 597623717, 528740935, 895670092, 741280370, 362867889, 849364135}
    };   
    }

    
    

    private void worstWafomsS3M5() {
    	System.out.println("for this methodes we construct the best one naively so we only take the result matrice as the best one");   }

    /*
     
     * M6:NAive sob+LMS 100 times
     * 
     * 
     * */
    private void meilleursWafomsS3M6() {
    	tabl = new int[][] {
    	{928245125,319456878,212602380,79922796,65572662,32679611,15717997,7732906,4013585,1990525,935470,447234,222241,118280,55749,24516,15090,7377,3879,1087},
    	{1058505031,770580776,856256822,604490373,1029770261,779778720,836252221,660563668,1059283383,767929827,859128065,606167146,1029125673,777150715,837170323,662359434,1058517184,770574518,856247146,604499227},
    	{988555712,421115755,791919253,818153746,393773032,677409388,1000894740,448188273,766325708,829199628,347994658,709286034,986193522,422346968,791497078,820426882,390148157,678768851,999192113,444789308}
    
    	};   }

    private void worstWafomsS3M6() {
    	tabl = new int[][] {
    		
    		{795480122,402101091,192074672,89722565,37018269,19937646,13535633,7447693,3437999,1738980,917070,458625,196019,114170,64592,29693,8196,7012,4033,1317},
    		{914229989,567284153,1055087242,770932271,897002993,589292948,1023920791,797411158,912167626,568737769,1053391522,768139509,894781792,589996870,1026523400,801089161,914243560,567295996,1055090175,770923871},
    		{1021552078,379676192,593781893,826722244,480879406,640501144,1036393587,350468605,545149099,815379050,512374141,631192670,1019818006,379235778,595162105,829078228,480107707,641254310,1039032429,349821013}
    		
    	};   }
    
    
    
    
    
    



    public int[][] getWafomsTabl() {
        return tabl;
    }

   /* public static void main(String[] args) {
        int s = 3;
        int methode=1;
        int BestOrWorstWaf=0;
        WafomsStorage storage = new WafomsStorage(s,methode,BestOrWorstWaf);
        int[][] meilleursWafoms = storage.getWafomsTabl();

        if (meilleursWafoms != null) {
            for (int i = 0; i < meilleursWafoms.length; i++) {
                for (int j = 0; j < meilleursWafoms[i].length; j++) {
                    System.out.print(meilleursWafoms[i][j] + " ");
                }
                System.out.println();
            }
        } else {
            System.out.println("ici on stock seulment les dimension s=3 et s=6 pour d'autres dimension faire la simulation sur ssj>src/main/docs/exemples/fichier_Sob.");
        }
    }*/
}