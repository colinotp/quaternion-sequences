if [ ! -d maplesat/maplesat ]
then
	git clone https://bitbucket.org/JLiangWaterloo/maplesat.git
fi
cd maplesat/maplesat
git checkout qw
make maplesat
cp simp/maplesat_static ../..
