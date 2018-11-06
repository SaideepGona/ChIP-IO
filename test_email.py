import smtplib

email = "ChIPBaseApp@gmail.com"
password = "chipbase"
target = "gona.saideep1@gmail.com"

def send_email(subject, msg):
    try:
        server = smtplib.SMTP('smtp.gmail.com:587')
        server.ehlo()
        server.starttls()
        server.login(email, password)
        message = 'Subject: {}\n\n{}'.format(subject, msg)
        server.sendmail(email, target, message)
        server.quit()
        print("Success: Email sent!")
    except:
        print("Email failed to send.")


subject = "Test subject"
msg = "Hello there, how are you today?"

send_email(subject, msg)