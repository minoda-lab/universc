version=1.2.7
old_version=1.2.6
sed -i "s/$old_version/$version/g"  *md *html .version man/*sh inst/* *.cff
sed -i 's/universcversion="$old_version"/universcversion="$version"/g' launch_universc.sh
R -e "knitr::knit('README.Rmd')"; pandoc -f markdown -t html README.md > README.html
echo "updating GitHub version $version"
git add -u
git reset HEAD test/
git commit -m "update documentation for $version"
git tag $version
git push --no-verify origin master
git push tag $version
echo "building Docker container version $version"
docker build -t universc:$version .
#docker tag  universc:$version universc:latest
#docker tag universc:$version docker.io/tomkellygenetics/universc:latest
docker tag universc:$version docker.io/tomkellygenetics/universc:$version
#docker push tomkellygenetics/universc:latest
docker push tomkellygenetics/universc:$version
exit 0
tag universc:$version docker.pkg.github.com/minoda-lab/universc/universc:$version
tag universc:$version docker.pkg.github.com/minoda-lab/universc/universc:latest
docker push docker.pkg.github.com/minoda-lab/universc/universc:$version
docker push docker.pkg.github.com/minoda-lab/universc/universc:latest
docker tag universc:$version docker.pkg.github.com/tomkellygenetics/universc/universc:$version
docker tag universc:$version docker.pkg.github.com/tomkellygenetics/universc/universc:latest
docker push docker.pkg.github.com/tomkellygenetics/universc/universc:latest
docker push docker.pkg.github.com/tomkellygenetics/universc/universc:$version
