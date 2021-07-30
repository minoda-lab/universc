version=1.1.4
old_version=1.1.3
sed -i "s/$old_version/$version/g"  *md *html launch_universc.sh .version man/*sh inst/*
R -e "knitr::knit('README.Rmd')"; pandoc -f markdown -t html README.md > README.html
echo "building version $version"
docker build -t universc:$version .
docker tag  universc:$version universc:latest
docker tag universc:$version tomkellygenetics/universc:latest
docker tag universc:$version tomkellygenetics/universc:$version
docker push tomkellygenetics/universc:latest
docker push tomkellygenetics/universc:$version
tag universc:$version docker.pkg.github.com/minoda-lab/universc/universc:$version
tag universc:$version docker.pkg.github.com/minoda-lab/universc/universc:latest
docker push docker.pkg.github.com/minoda-lab/universc/universc:$version
docker push docker.pkg.github.com/minoda-lab/universc/universc:latest
docker tag universc:$version docker.pkg.github.com/tomkellygenetics/universc/universc:$version
docker tag universc:$version docker.pkg.github.com/tomkellygenetics/universc/universc:latest
docker push docker.pkg.github.com/tomkellygenetics/universc/universc:latest
docker push docker.pkg.github.com/tomkellygenetics/universc/universc:$version
