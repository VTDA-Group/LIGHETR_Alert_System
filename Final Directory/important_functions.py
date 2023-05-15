# https://www.youtube.com/watch?v=g_j6ILT-X0k

from email.message import EmailMessage
import ssl
import smtplib
import mimetypes
import healpy as hp



def dial_numbers(numbers_list):
    """Dials one or more phone numbers from a Twilio phone number."""
    for number in numbers_list:
        print("Dialing " + number)
        # set the method to "GET" from default POST because Amazon S3 only
        # serves GET requests on files. Typically POST would be used for apps
        client.calls.create(to=number, from_=TWILIO_PHONE_NUMBER,
                            url=TWIML_INSTRUCTIONS_URL, method="GET")



def send_email(email_sender, email_password, all_email_recipients, files, subject, body):
    
    print("Sending email, sender: "+str(email_sender)+" pass: "+str(email_password)+" recipients: "+str(all_email_recipients))
    
    everyone = " , ".join(all_email_recipients)
    subject = subject
    body = body+"\nList of all Recipients of this email: "+str(everyone)
    
    for recipients in all_email_recipients:
        
        print("recipients: "+str(recipients))
        
        
        em = EmailMessage()
        em['From'] = email_sender
        em['To'] = recipients
        em['Subject'] = subject
        em.set_content(body)

        
        #attaching files
        for path in files:
            ctype, encoding = mimetypes.guess_type(path)
            if ctype is None or encoding is not None:
                # No guess could be made, or the file is encoded (compressed), so
                # use a generic bag-of-bits type.
                ctype = 'application/octet-stream'
            maintype, subtype = ctype.split('/', 1)

            with open(path, 'rb') as fp:
                em.add_attachment(fp.read(), maintype = maintype, subtype = subtype, filename = path)

        '''
        #sending email
        context = ssl.create_default_context()

        with smtplib.SMTP_SSL("smtp.gmail.com", 465, context = context) as smtp:
            smtp.login(email_sender, email_password)
            smtp.sendmail(email_sender, recipients, em.as_string())'''
