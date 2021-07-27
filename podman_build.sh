version=1.1.4
sed -i "s/1.1.3/1.1.4/g"  *md *html  launch_universc.sh .version man/*sh inst/*
R -e "knitr::knit('README.Rmd')"; pandoc -f markdown -t html README.md > README.html
echo "building version $version"
podman build -t universc:1.1.3 .
podman tag  universc:1.1.3 universc:latest
podman tag universc:1.1.3 tomkellygenetics/universc:latest
podman tag universc:1.1.3 tomkellygenetics/universc:1.1.3
podman push tomkellygenetics/universc:latest
podman push tomkellygenetics/universc:1.1.3
podman tag universc:1.1.3 docker.pkg.github.com/minoda-lab/universc/universc:1.1.3
podman tag universc:1.1.3 docker.pkg.github.com/minoda-lab/universc/universc:latest
podman push docker.pkg.github.com/minoda-lab/universc/universc:latest
podman push docker.pkg.github.com/minoda-lab/universc/universc:1.1.3
podman tag universc:1.1.3 docker.pkg.github.com/tomkellygenetics/universc/universc:1.1.3
podman tag universc:1.1.3 docker.pkg.github.com/tomkellygenetics/universc/universc:latest
podman push docker.pkg.github.com/tomkellygenetics/universc/universc:latest
podman push docker.pkg.github.com/tomkellygenetics/universc/universc:1.1.3
