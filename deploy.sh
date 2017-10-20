./build-doc
echo
(
    cd ..
    rm -f AnchoredInversionClient_*.tar.gz
    echo
    R CMD REMOVE AnchoredInversionClient
    echo
    R CMD build client.R
    echo
    R CMD INSTALL AnchoredInversionClient_*.tar.gz
)

