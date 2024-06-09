### Streamlit
from random import shuffle

import streamlit as st
import streamlit.components.v1 as components

from methyl_aware_distance import MethylAwareDistance


import pandas as pd
import numpy as np
import altair as alt


#### Styling functions
stHTML = lambda string: st.markdown(string, unsafe_allow_html=True)
stHTMLsidebar = lambda string: st.sidebar.markdown(string, unsafe_allow_html=True)
textCenter = lambda text: f"<p style='text-align:center;'> {text} </p>"
center_string = lambda x: f"<p style='text-align:center;'> {x} </p>"


##########################################

########################################## SIDEBAR #############################################
#st.sidebar.title('Sidebar')
st.sidebar.markdown('---')


html_string = "<h1 style='text-align:center;'> An interactive tool for methyl-aware edit distance calculations </h1>"
stHTML(html_string)
st.markdown('---')




expected_sequence = st.text_input("Expected sequence","ACTG")
observed_sequence = st.text_input("Observed sequence","ATTG")



# Boolean input (checkbox)
is_read_1 = st.checkbox("Sequence from Read 1? (if read 2 select False)",True)


if expected_sequence and observed_sequence:

    edit_distance = MethylAwareDistance(expected_sequence,
                                        observed_sequence,
                                        is_read_1)
    
    if edit_distance.equal_length:
        st.write("Hamming distance", edit_distance.hamming_methyl_aware())
    else:
        st.write("Levenshtein distance", edit_distance.space_efficient_levenshtein_methyl_aware())



# #st.sidebar.markdown("<h1 style='text-align:center;'> This title is HTML </h1>", unsafe_allow_html=True)

# stHTMLsidebar("<h1 style='text-align:center;'> Distribution parameters </h1>")

# ## Selector
# st.sidebar.markdown('Select average weight :rabbit:  :point_down:')
# mu = st.sidebar.slider('',1,6,value=3,step=1)

# stHTMLsidebar('<br>')

# sigma = st.sidebar.slider('+-:',0.0,3.0,value=0.5,step=0.5)


# st.markdown('<br>', unsafe_allow_html=True)
# stHTML('<br>')

# stHTMLsidebar("<h1 style='text-align:center;'> Sample size </h1>")

# st.sidebar.markdown('Select sample size of :rabbit2: weights :point_down:')
# sample_size = st.sidebar.slider('sample size:',1,30,value=3,step=1)

# st.sidebar.markdown('Number of random samples :point_down:')
# n_samples = st.sidebar.slider('sample size:',1,100,value=5,step=1)

# st.sidebar.markdown('---')

# #############################################################################################


# ########################################## MAIN #############################################

# html_string = "<h1 style='text-align:center;'> Sampler with normal distribution </h1>"
# stHTML(html_string)
# st.markdown('---')



# #############################################


# def get_pdf(mu, sigma):
#     # Mean and standard deviation

#     # Generating points on the x axis between -4 and 4, which is typical for a normal distribution plot
#     x = np.linspace(mu - 4*sigma, mu + 4*sigma, 1000)

#     # Calculating the probability density function for the normal distribution
#     pdf = (1 / (np.sqrt(2 * np.pi * sigma**2))) * np.exp(-0.5 * ((x - mu) ** 2) / (sigma ** 2))
#     return pd.DataFrame({'X':x,'PDF':pdf})


# # Step 2: Function to extract a sample of n weights from the distribution
# def sample_weights(mu, sigma, n):
#     return np.random.normal(mu, sigma, n)

# # Step 3: Function to plot the mean of the sample onto the previous plot
# def get_sample_mean(mu, sigma, n):
#     sample = sample_weights(mu, sigma, n)
#     return np.mean(sample)


# underlying_distribution = get_pdf(mu, sigma)


# sample_means = [(get_sample_mean(mu,sigma, sample_size),.01) for i in range(n_samples)]
# sample_means = pd.DataFrame(sample_means,columns=['mean','ypos'])


# #Chart with tooltip 
# st.markdown('<style>#vg-tooltip-element{z-index: 1000051}</style>', unsafe_allow_html=True)#hack to make tooltip visible in full-screen mode

# chart = alt.Chart(underlying_distribution).mark_line().encode(
#     x=alt.X('X',title=f'distribution of weight'),
#     y=alt.Y('PDF',axis=None),
#     color=alt.value('black')

# ).properties(
#     title=f'Normal Distribution of Rabbit Weights | Mean={mu} kg, SD={sigma} kg',
#     width=600,
#     height=300,
    
# )



# points_chart = alt.Chart(sample_means).mark_point(filled=True).encode(
#     x='mean',
#     y='ypos',
#     color=alt.value('red')
# )


# chart = (chart + points_chart).interactive()
# chart = chart.configure(background='white')
# st.altair_chart(chart, use_container_width=True)



# ### Some space
# stHTML('<br>')


#############################################
stHTML('<br>')
stHTML('<br>')
st.markdown('---')
st.write('Joel Rodriguez-Medina')
#logoDir = 'https://1jj36121ulnl3yz5ue21q04x-wpengine.netdna-ssl.com/wp-content/uploads/2020/05/Foodome_Logo.png'
#components.html(f"<img src='{logoDir}' style='width:25%;margin-left: auto;margin-right: auto;display: block'>")