#!/usr/bin/env python
# coding: utf-8

# <div align='center' style='font-size:100%'><b>Biological networks and data aggregation</b></div>
# <br>
# <div align='center'><img src="https://raw.githubusercontent.com/vporubsky/tellurium-libroadrunner-tutorial/master/data_aggregation_logo.png" width = "50%" style="padding: 0px"></div>
# <br>
# <div align='center' style='font-size:100%'>
# Veronica L. Porubsky, BS
# <div align='center' style='font-size:100%'>Sauro Lab PhD Student, Department of Bioengineering<br>
# Head of Outreach, <a href="https://reproduciblebiomodels.org/dissemination-and-training/seminar/">Center for Reproducible Biomedical Modeling</a><br>
# University of Washington, Seattle, WA USA
# </div>
# <hr>

# ## TOC
# * [What is biochemical network modeling](#network-modeling)
# * [Why do we perform network modeling](#perform-network-modeling)
# * [Types of networks](#network-types)
# * [Metabolic networks](#metabolic)
# * [Protein signaling networks](#protein-signaling)
# * [Gene regulatory networks](#gene-regulation)
# * [Repressilator model by Elowitz & Liebler (2000):](#repressilator)
# * [Model of respiratory oscillations in Saccharomyces cerevisae by Jana Wolf et al. (2001):](#wolf)
# * [Negative feedback and ultrasensitivity effects in Kholodenko (2000):](#kholodenko)
# * [What is data aggregation?](#data-aggregation)
# * [What databases are useful for biochemical network modeling?](#databases)
# * [What is metadata and how much should we collect?](#metadata)
# * [Packages and Constants](#packages-constants)
# * [Importing data with KEGG and bioservices](#import-with-KEGG)
# * [Aggregating data with ChEBI](#aggregation-with-chebi)
# * [Exercises](#exercises)

# #  What is biochemical network modeling? <a class="anchor" id="network-modeling"></a>
# 
# <ul>
#   <li>Chemical kinetics studies the factors that influence the rate of chemical reactions</li>
#      <ul class="square">
#       <li>e.g. <span style="color:blue">concentration</span>, temperature, light, catalysts, etc. </li>
#      </ul>
#   <li>Chemical reaction networks are the framework for building all types of dynamical models</li>
#         <ul class="square">
#           <li>Genetic circuits</li>
#           <li>Cell signaling pathways</li>
#           <li>Metabolic networks</li>
#         </ul>
#   <li>Types of biochemical network models:</li>
#     <ul class="square">
#       <li>Agent-based</li>
#       <li>Algebraic</li>
#       <li>Boolean</li>
#       <li>Constraint based</li>
#         <li><span style="color:blue">Mechanistic differential equations models</span></li>
#       <li>Statistical and machine learning methods</li>
#       <li>Stochastic</li>
#     </ul>
# </ul>

# #  Why do we perform network modeling? <a class="anchor" id="perform-network-modeling"></a>
# 
# <ul>
#   <li>Understand subcellular processes</li>
#   <li>Drive experimentation</li>
#   <li>Make predictions about system behavior and the impacts of interferring with the system</li>
#   <li>Design synthetic networks to control cellular processes</li>
#   <li>Develop novel treatments for disease by predicting targets for pharmacological therapies, genetic modification, etc</li>
#   <li>Provide a basis for larger multi-cellular models</li>
# </ul>
# 

# # Types of networks: <a class="anchor" id="network-types"></a>
# 
# <br>
# <div align='center'><img src="https://raw.githubusercontent.com/vporubsky/tellurium-libroadrunner-tutorial/master/networks_fig1.PNG" width = "75%" style="padding: 0px"></div>
# <br>

# # Metabolic networks: <a class="anchor" id="metabolic"></a>
# 
# <br>
# Glycolytic pathway of <em>Lactococcus lactis</em>, produced on JWS Online for Dr. Sauro's textbook, "Systems Biology: Introduction to Pathway Modeling".
# <br>
# 
# <table><tr>
#    <td> <img src="https://raw.githubusercontent.com/vporubsky/tellurium-libroadrunner-tutorial/master/metabolic_network.PNG" style="width: 800px;"/> </td>
# </tr></table>
# 
# 

# # Protein signaling networks: <a class="anchor" id="protein-signaling"></a>
# 
# <table><tr>
#    <td> <img src="https://raw.githubusercontent.com/vporubsky/tellurium-libroadrunner-tutorial/master/protein_actions.PNG" style="width: 450px;"/> </td>
#    <td> <img src="https://raw.githubusercontent.com/vporubsky/tellurium-libroadrunner-tutorial/master/protein_signalling_network.PNG" style="width: 300px;"/> </td>
# </tr></table>

# # Gene regulatory networks: <a class="anchor" id="gene-regulation"></a>
# 
# <table><tr>
#    <td> <img src="https://raw.githubusercontent.com/vporubsky/tellurium-libroadrunner-tutorial/master/gene_actions.PNG" style="width: 400px;"/> </td>
#    <td> <img src="https://raw.githubusercontent.com/vporubsky/tellurium-libroadrunner-tutorial/master/gene_regulatory.PNG" style="width: 400px;"/> </td>
# </tr></table>

# #  Repressilator model by Elowitz & Liebler (2000): <a class="anchor" id="repressilator"></a>
# 
# <br>
# Repressilator circuit from <a href="http://www.elowitz.caltech.edu/publications/Repressilator.pdf">Elowitz & Liebler (2000):</a>
# <div align='center'><img src="https://raw.githubusercontent.com/vporubsky/tellurium-libroadrunner-tutorial/master/repressilator.png" width="50%" style="padding: 20px"></div>
# 

# #  Repressilator model by Elowitz & Liebler (2000):
# 
# <table><tr>
#    <td> <img src="https://raw.githubusercontent.com/vporubsky/tellurium-libroadrunner-tutorial/master/repressilator_2.PNG" style="width: 300px;"/> </td>
#    <td> <img src="https://raw.githubusercontent.com/vporubsky/tellurium-libroadrunner-tutorial/master/repressilator_3.PNG" style="width: 1000px;"/> </td>
# </tr></table>

# #  Repressilator model by Elowitz & Liebler (2000):
# 
# <table><tr>
#    <td> <img src="https://raw.githubusercontent.com/vporubsky/tellurium-libroadrunner-tutorial/master/repressilator_1.PNG" style="width: 700px;"/> </td>
# </tr></table>

# #  Repressilator model by Elowitz & Liebler (2000):
# 
# <table><tr>
#    <td> <img src="https://raw.githubusercontent.com/vporubsky/tellurium-libroadrunner-tutorial/master/repressilator_4.PNG" style="width: 700px;"/> </td>
# </tr></table>

# # Model of respiratory oscillations in Saccharomyces cerevisae by Jana Wolf et al. (2001): <a class="anchor" id="wolf"></a>
# 
# <div align='center'><img src="https://raw.githubusercontent.com/vporubsky/tellurium-libroadrunner-tutorial/master/wolf_publication.PNG" width="65%" style="padding: 20px"></div>
# <br>
# <div align='center'><img src="https://raw.githubusercontent.com/vporubsky/tellurium-libroadrunner-tutorial/master/wolf_network.PNG" width="65%" style="padding: 20px"></div>

# # Negative feedback and ultrasensitivity effects in Kholodenko (2000): <a class="anchor" id="kholodenko"></a>
# 
# <table><tr>
#    <td> <img src="https://raw.githubusercontent.com/vporubsky/tellurium-libroadrunner-tutorial/master/kholodenko_1.PNG" style="width: 500px;"/> </td>
#    <td> <img src="https://raw.githubusercontent.com/vporubsky/tellurium-libroadrunner-tutorial/master/kholodenko_2.PNG" style="width: 500px;"/> </td>
# </tr></table>
# 
# 

# # What is data aggregation? <a class="anchor" id="data-aggregation"></a>
# 
# The collection of data from multiple experiments, scientific papers, and online data sources. Typically you will need to curate your aggregated data to ensure the quality of measurements are suffiencient to include in your model.
# 

# # What databases are useful for biochemical network modeling? <a class="anchor" id="databases"></a>
# 
# 
# <ul>
#   <li>SABIO-RK: biochemical reaction kinetics database</li>
#      <ul class="square">
#       <li>Describes chemical reactions and kinetics</li>
#       <li>Contains information about participants and modifiers in reactions</li>
#       <li>Metabolic and signaling network reactions</li>
#      </ul>
#   <li>BRENDA: the comprehensive enzyme information system</li>
#      <ul class="square">
#       <li>Enzyme information classified by the biochemical reaction it catalyzes</li>
#       <li>Kinetic information about substrates and products is available</li> 
#      </ul>
#   <li>ChEBI: dictionary of "small" chemical compounds</li>
#   <li>KEGG: collection of pathway/genome/diesease/drug databases</li>
#   <li>BioCYC: collection of pathway/genome databases</li>
#        <ul class="square">
#       <li>Search for genes, proteins, metabolites or pathways, and the occurence of your term will be located in multiple databases</li> 
#      </ul>
#   <li>BioModel: repository of mathematical models of biological systems</li>
#       <ul class="square">
#       <li> *Will be covered in more detail later in the course</li> 
#      </ul>
# </ul>
# 
# <a href="https://www.sciencedirect.com/science/article/abs/pii/S0958166917301428?via%3Dihub">Appendix A of Goldberg et al. (2018)</a> provides a useful and more comprehensive list of data sources containing intracellular biochemical data. 
# 

# # What is metadata and how much should we collect? <a class="anchor" id="metadata"></a>
# 
# <ul>
#   <li>Metadata: data that describes biochemical data</li>
#   <li>Collect information about:</li>
#      <ul class="square">
#       <li>Units</li>
#       <li>Estimates of measurement accuracy</li>
#       <li>Annotations</li>
#       <li>Ontology terms defining the annotations</li>
#       <li>etc.</li>
#      </ul>
#   <li>Collect provenance data:</li>
#      <ul class="square">
#       <li>Lab which generated the data</li>
#       <li>Experimental conditions</li>
#       <li>Protocol used to generate the data</li>
#       <li>Paper which reported the measurement</li>
#       <li>etc.</li>
# </ul>

# # Packages and constants <a class="anchor" id="packages-constants"></a>

# In[9]:


from bioservices import *


# # Importing data with KEGG and bioservices <a class="anchor" id="import-with-KEGG"></a>

# In[10]:


# Select database
database = KEGG()

# Retrieve a KEGG entry
tetR_query = database.get("K18476")

# Build a dictionary to parse query
tetR_dict = database.parse(tetR_query)

# Show information about the query
print(tetR_dict['NAME'])
print(tetR_dict['DEFINITION'])


# # Aggregating data with ChEBI <a class="anchor" id="aggregation-with-chebi"></a>

# In[11]:


# Store annotation information
# Select database
database = ChEBI()

# Retrieve a ChEBI entry for D-fructose 1,6-bisphosphate
query = database.getCompleteEntity("CHEBI:78682")

print(query.definition)


# # Exercises <a class="anchor" id="exercises"></a>

# ## Exercise 1
# 
# Visit SABIO-RK and find a reaction involving tau-protein. What tissue and organism is the provided reaction relevant to, according to the results?
# 
# **Solution:** Tissue: brain, organisms: Homo sapiens and Rattus norvegicus 
# 
# What other metadata are available for the reaction
# 
# **Solution:**
# 

# ## Exercise 2
# 
# Visit BRENDA and search for the enzyme lactase. What reaction does this enzyme catalyze? Record the reactants and products in the reaction.
# 
# 
# **Solution:**
# 
# What other databases does BRENDA link to which could help you build a pathway model containing the specified reaction?
# 
# 
# **Solution:**

# ## Exercise 3
# 
# Craving a coffee? Look up "caffeine" in ChEBI. Getting late in your region of the world? Try "melatonin" instead.
# 
# **Solution:**
# 
# What is the ChEBI ID for your small molecule? What organisms was the metabolite detected in or isolated from according to the ChEBI data?
# 
# **Solution:**
# 

# ## Exercise 4
# 
# Explore this <a href="https://biocyc.org/overviewsWeb/celOv.shtml?orgid=ECOLI"> full metabolic map on BioCYC. </a> 

# ## Exercise 5
# 
# Programmatically access the ChEBI entry for glucose and print the molecular formula information.

# In[12]:


# Exercise 5 Solution:

# Select database
database = ChEBI()

# Retrieve a ChEBI entry for D-fructose 1,6-bisphosphate
query = database.getCompleteEntity("CHEBI:17234")

print(query.Formulae)


# # Acknowledgements
# <br>
# <div align='left'><img src="https://raw.githubusercontent.com/vporubsky/tellurium-libroadrunner-tutorial/master/acknowledgments.png" width="80%"></div>

# <br>
# <html>
#    <head>
#       <title>Bibliography</title>
#    </head>
#    <body>
#       <h1>Bibliography</h1>
#       <ol>
#          <li>
#             <p>K. Choi et al., <cite>Tellurium: An extensible python-based modeling environment for systems and synthetic biology</cite>, Biosystems, vol. 171, pp. 74–79, Sep. 2018.</p>
#          </li>
#          <li>
#             <p>E. T. Somogyi et al., <cite>libRoadRunner: a high performance SBML simulation and analysis library.,</cite>, Bioinformatics, vol. 31, no. 20, pp. 3315–21, Oct. 2015.</p>         
#           <li>
#             <p>L. P. Smith, F. T. Bergmann, D. Chandran, and H. M. Sauro, <cite>Antimony: a modular model definition language</cite>, Bioinformatics, vol. 25, no. 18, pp. 2452–2454, Sep. 2009.</p>
#          </li>
#          <li>
#             <p>K. Choi, L. P. Smith, J. K. Medley, and H. M. Sauro, <cite>phraSED-ML: a paraphrased, human-readable adaptation of SED-ML</cite>, J. Bioinform. Comput. Biol., vol. 14, no. 06, Dec. 2016.</p>
#          </li>         
#          <li>
#             <p> B.N. Kholodenko, O.V. Demin, G. Moehren, J.B. Hoek, <cite>Quantification of short term signaling by the epidermal growth factor receptor.</cite>, J Biol Chem., vol. 274, no. 42, Oct. 1999.</p>
#          </li>
#       </ol>
#    </body>
# </html>
