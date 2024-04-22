# Ecological Network Inference Validation
Part of Chapter 3 of my PhD project.

## Abstract

The labour-intensive sampling requirements of ecological networks have spurred the creation of network inference methodology. Demonstrated inconsistencies in networks inference necessitates quantification of inference performance to facilitate choice of network inference approaches and parameterisation. Here, we present a simulation-validation framework which generates data products fit for network inference and subsequently assesses inference performance. The simulation framework we present here can be parameterised using real-world data. (e.g., biological interactions observed in-situ and bioclimatic niche preferences). In addition, our simulation framework supports directed and undirected links in ecological networks as well as exporting of time-series or spatial products fit for currently available ecological network inference approaches. The validation procedures we present quantify inference and detection probabilities of association types of different identity and sign with respect to bioclimatic niche preferences and association strength between species. Applying our workflow to one well-established ecological association network inference method (HMSC), we identify a large range in accuracy of inferred networks as compared to true, realisable association networks. We find that these differences in inference accuracy are governed by a paradigm of input data types and environmental parameter estimation with network inference approach parameterisation accounting for environmental gradients and leveraging more nuanced biological inputs outperforming simpler parameterisations in inference accuracy. Conclusively, we provide the groundwork for validation and comparison of ecological network inference approaches to improve our capabilities of predicting species biodiversity across space and time.

## Research Questions
How performant are network inference approaches and what affects their accuracy?

## Primary Contact and Collaborators
### Primary Contact
Erik Kusch (erik.kusch@nhm.uio.no, [OrcID](https://orcid.org/my-orcid?orcid=0000-0002-4984-7646))  

### Collaborators
- Anna C. Vinton ([OrcID](https://orcid.org/0000-0002-8279-1736))

## Key Findings
1. Using co-performance or co-abundance rather than the traditional co-occurrence inputs improves inference accuracy
2. Accounting for environmental gradients improves inference accuracy

Network inference approaches have been developed to address the labour intensity of ecological network research which requires impractical sampling efforts particularly at large geographical scales or for large species pools. Despite the large number and diversity in approaches of network inference frameworks, no general framework for network inference validation or performance quantification has been established yet. Here, we have presented an easy-to-use, mathematically robust simulation framework which can be used for data generation with which to carry out network inference validation. While we find promising accuracy of network inference under specific conditions (i.e., use of site-species performance matrices and environmental information), we argue that these performance scores may be subject to change given further parametrising of our simulation with more realistic real-world data. Additionally, further explorations of network inference accuracy of directed links (i.e., interactions) or using the capacity of our simulation framework to output time-series of individuals across space may prove particularly useful in aiding to the overview of network inference capabilities and current practices thereof.
Nevertheless, using our simulation framework and a novel network inference performance and behaviour quantification pipeline, we have confirmed previous criticisms of ecological network inference practices and added to existing guidelines for choice of appropriate network inference approaches. Thus, our simulation framework for network inference validation represents an important addition to the toolbox of network inference method development which fills a crucial knowledge gap of network inference performance assessments.

## Data Availability
Data needed for reproducing this analysis is available at https://zenodo.org/doi/10.5281/zenodo.10535149.

## Funding
Aarhus University Research Foundation Start-up Grant (grant no. AUFF-2018-7-8) 
European Union Horizon Europe (grant no. [101057437](https://doi.org/10.3030/101057437))
