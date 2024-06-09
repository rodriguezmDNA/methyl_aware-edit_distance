### Streamlit

import streamlit as st

from methyl_aware_distance import MethylAwareDistance


#### Styling functions
stHTML = lambda string: st.markdown(string, unsafe_allow_html=True)
stHTMLsidebar = lambda string: st.sidebar.markdown(string, unsafe_allow_html=True)
textCenter = lambda text: f"<p style='text-align:center;'> {text} </p>"
center_string = lambda x: f"<p style='text-align:center;'> {x} </p>"


##########################################

########################################## SIDEBAR ##########################################
#st.sidebar.title('Sidebar')
#st.sidebar.markdown('---')


###########################################  MAIN  ###########################################
html_string = "<h2 style='text-align:center;'> A tool for methyl-aware edit distance calculations </h2>"
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

        html_string = f"<h4 style='text-align:center;'> Hamming distance: <code>{edit_distance.hamming_methyl_aware()}</code> </h4>"
        

        #st.write("Hamming distance", edit_distance.hamming_methyl_aware())
    else:
        html_string = f"<h4 style='text-align:center;'> Levenshtein distance: <code>{edit_distance.space_efficient_levenshtein_methyl_aware()}</code></h4>"
        #st.write("Levenshtein distance", edit_distance.space_efficient_levenshtein_methyl_aware())
    st.write('\n')
    st.write('\n')
    st.write('\n')
    stHTML(html_string)


#############################################
stHTML('<br>')
stHTML('<br>')
st.markdown('---')


html_string = f"<h5 style='text-align:center;'> Joel Rodriguez Medina </h5>"
stHTML(html_string)