#include "testkeycomponents.h"

testkeycomponents::testkeycomponents(size_t k_, percent_t pkps)
{
	k = k_;
	percent_kps = pkps;
}

void testkeycomponents::configure(size_t k_, percent_t pkps)
{
	k = k_;
	percent_kps = pkps;
}

testkeycomponents::~testkeycomponents()
{

}

void testkeycomponents::run(string database, string sub)
{
	string file;
	ifstream is(PATH_DATA + database);
	while(is>>file)
		one_test_fm(file, sub); //test type

	is.close();
}
void testkeycomponents::one_test(string file, string sub)
{
	off shape(PATH_DATA + file);
	size_t nkps = shape.get_nvertices()*percent_kps;

	keypoints kps(shape);

	keycomponents kcs(kps,shape,k,percent_kps);

	ofstream oskc(PATH_KCS + sub + "/"+ getFileName(file) + "kc");
	kcs.print(oskc);
	oskc.close();

	//ofstream oskp(PATH_KPS + sub + "/"+ getFileName(file) + "kp");
	//kps.print(oskp, nkps);
	//oskp.close();

	//NEW KPS
	ofstream oskp(PATH_KPS + sub + "/"+ getFileName(file) + "kp");
	kcs.new_keypoints(shape, kps, oskp, 0.1);
	oskp.close();

	cout<<"INCUNTRU: "<<kcs.get_ncomponents()<<endl;
}

void testkeycomponents::one_test_fm(string file, string sub)
{
	cout<<"KPS"<<endl;
	off shape(PATH_DATA + file);

	size_t nkps = shape.get_nvertices()*percent_kps;

	keypoints kps(shape);
	ofstream oskp(PATH_KPS + sub + "/"+ getFileName(file) + "kp");
	kps.print(oskp, nkps);
	oskp.close();


	vector<index_t> source(kps.get_keypoints(), kps.get_keypoints() + nkps);
	fastmarching fm(shape, source, INFINITY, 1);

	ofstream os(PATH_TEST + string("test.fm"));
	fm.print(os);
	os.close();

	percent_t radio = k;
	radio /= 100;

	keycomponents kcs(kps, shape, fm, radio, percent_kps);

	ofstream oskc(PATH_KCS + sub + "/"+ getFileName(file) + "kc");
	kcs.print_fm(oskc);
	oskc.close();
	cout<<"INCUNTRU: "<<kcs.get_ncomponents()<<endl;
}
