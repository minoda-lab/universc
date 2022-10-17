version=1.2.4
old_version=1.2.3
sed -i "s/$old_version/$version/g"  *md *html  launch_universc.sh .version man/*sh inst/* *.cff
#R -e "knitr::knit('README.Rmd')"; pandoc -f markdown -t html README.md > README.html
echo "updating GitHub version $version"
git add -u
git reset HEAD test/
git commit -m "update documentation for $version"
git tag $version
git push --no-verify origin master
git push tag $version
echo "building Docker container version $version"
podman build -t universc:$version .
podman tag  universc:$version universc:latest
podman tag universc:$version docker.io/tomkellygenetics/universc:latest
podman tag universc:$version docker.io/tomkellygenetics/universc:$version
podman push docker.io/tomkellygenetics/universc:latest
podman push docker.io/tomkellygenetics/universc:$version
podman tag universc:$version docker.pkg.github.com/minoda-lab/universc/universc:$version
podman tag universc:$version docker.pkg.github.com/minoda-lab/universc/universc:latest
podman push docker.pkg.github.com/minoda-lab/universc/universc:latest
podman push docker.pkg.github.com/minoda-lab/universc/universc:$version
podman tag universc:$version docker.pkg.github.com/tomkellygenetics/universc/universc:$version
podman tag universc:$version docker.pkg.github.com/tomkellygenetics/universc/universc:latest
podman push docker.pkg.github.com/tomkellygenetics/universc/universc:latest
podman push docker.pkg.github.com/tomkellygenetics/universc/universc:$version
